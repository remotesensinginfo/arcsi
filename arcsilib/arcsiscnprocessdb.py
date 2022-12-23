import logging
import os
import shutil
from sqlite3 import Connection as SQLite3Connection

import sqlalchemy
import sqlalchemy.orm
from sqlalchemy import event
from sqlalchemy.engine import Engine
from sqlalchemy.ext.declarative import declarative_base

logger = logging.getLogger(__name__)

Base = declarative_base()


class ARCSIScnProcess(Base):
    __tablename__ = "ARCSIScnProcess"
    product_id = sqlalchemy.Column(sqlalchemy.String, primary_key=True)
    sensor = sqlalchemy.Column(sqlalchemy.String, primary_key=True)
    scn_url = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    geo_str_id = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    download = sqlalchemy.Column(sqlalchemy.Boolean, nullable=False, default=False)
    download_path = sqlalchemy.Column(sqlalchemy.String, nullable=True)
    ard = sqlalchemy.Column(sqlalchemy.Boolean, nullable=False, default=False)
    ard_path = sqlalchemy.Column(sqlalchemy.String, nullable=True)


@event.listens_for(Engine, "connect")
def set_sqlite_pragma(dbapi_connection, connection_record):
    if isinstance(dbapi_connection, SQLite3Connection):
        cursor = dbapi_connection.cursor()
        cursor.execute("PRAGMA journal_mode = MEMORY")
        cursor.execute("PRAGMA synchronous = OFF")
        cursor.execute("PRAGMA temp_store = MEMORY")
        cursor.execute("PRAGMA cache_size = 500000")
        cursor.close()


class RecordScn2Process(object):
    def __init__(self, sqlite_db_file):
        """
        Constructor for the class.

        :param sqlite_db_file: A file path for the SQLite database must be provided.

        """
        self.sqlite_db_file = sqlite_db_file
        self.sqlite_db_conn = "sqlite:///{}".format(self.sqlite_db_file)

    def init_db(self):
        """
        A function which must be called before use if a database file does not already
        exist. Note. if the database does exist then it will be deleted and recreated
         and any data with the existing database will be lost.

        """
        try:
            logger.debug("Creating Database Engine and Session.")
            db_engine = sqlalchemy.create_engine(
                self.sqlite_db_conn, pool_pre_ping=True
            )
            Base.metadata.drop_all(db_engine)
            logger.debug("Creating Database.")
            Base.metadata.bind = db_engine
            Base.metadata.create_all()
            session_sqlalc = sqlalchemy.orm.sessionmaker(bind=db_engine)
            ses = session_sqlalc()
            ses.close()
            logger.debug("Created Database Engine and Session.")
        except:
            raise Exception(
                "The SQLite database file cannot be opened: '{}'".format(
                    self.sqlite_db_conn
                )
            )

    def add_scns(self, scns_lst):
        """
        A function which adds scenes to the database.

        The geo_str_id can largely be what ever you want to use but for landsat, this is probably 'r{row}_p{path}'
        and for sentinel-2 then it is the granule.

        :param scns_lst: a list of dicts where the input dict must contain the following keys:
                         'product_id', 'sensor', 'scn_url', 'geo_str_id'.

        """
        logger.debug("Creating Database Engine.")
        db_engine = sqlalchemy.create_engine(self.sqlite_db_conn, pool_pre_ping=True)
        logger.debug("Creating Database Session.")
        session_sqlalc = sqlalchemy.orm.sessionmaker(bind=db_engine)
        ses = session_sqlalc()
        logger.debug("Created Database Engine and Session.")

        scn_lst = []
        for scn in scns_lst:
            scn_lst.append(
                ARCSIScnProcess(
                    product_id=scn["product_id"],
                    sensor=scn["sensor"],
                    scn_url=scn["scn_url"],
                    geo_str_id=scn["geo_str_id"],
                )
            )

        logger.debug(
            "There are {} scenes to be written to the database.".format(len(scn_lst))
        )
        if len(scn_lst) > 0:
            ses.add_all(scn_lst)
            ses.commit()
            logger.debug("Written jobs to the database.")
        ses.close()

    def is_scn_in_db(self, product_id, sensor):
        """
        A function to test whether a scene is within the database.

        :param product_id: the product ID for the scene
        :param sensor: the sensor for the scene
        :return: boolean (True: is present within the database. False: is not presented in the database).

        """
        logger.debug("Creating Database Engine.")
        db_engine = sqlalchemy.create_engine(self.sqlite_db_conn, pool_pre_ping=True)
        logger.debug("Creating Database Session.")
        session_sqlalc = sqlalchemy.orm.sessionmaker(bind=db_engine)
        ses = session_sqlalc()
        logger.debug("Created Database Engine and Session.")

        logger.debug("Perform query to find scene.")
        query_result = (
            ses.query(ARCSIScnProcess)
            .filter(
                ARCSIScnProcess.product_id == product_id,
                ARCSIScnProcess.sensor == sensor,
            )
            .one_or_none()
        )
        ses.close()
        logger.debug("Closed the database session.")
        found_prod_id = True
        if query_result is None:
            found_prod_id = False
        return found_prod_id

    def n_geoid_scns(self, geo_str_id):
        """
        A function to get a count of the number of scenes available in the database for a geographic filter.

        :param geo_str_id: A geographic filter which is dependent on the sensor.
                           If landsat, this is 'r{row}_p{path}. If sentinel-2 then it is the granule.'
                           Default: None (i.e., no filter applied).
        :return: int

        """
        logger.debug("Creating Database Engine.")
        db_engine = sqlalchemy.create_engine(self.sqlite_db_conn, pool_pre_ping=True)
        logger.debug("Creating Database Session.")
        session_sqlalc = sqlalchemy.orm.sessionmaker(bind=db_engine)
        ses = session_sqlalc()
        logger.debug("Created Database Engine and Session.")

        logger.debug("Perform query to find scene.")
        n_scns = (
            ses.query(ARCSIScnProcess)
            .filter(ARCSIScnProcess.geo_str_id == geo_str_id)
            .count()
        )
        ses.close()
        logger.debug("Closed the database session.")
        return n_scns

    def geoid_scns(self, geo_str_id):
        """
        A function which gets a list of scenes for a geographic filter (dependent on sensor). This function
        returns a list of scenes regardless of the processing stage.

        :param geo_str_id: A geographic filter which is dependent on the sensor.
                           If landsat, this is 'r{row}_p{path}. If sentinel-2 then it is the granule.'
                           Default: None (i.e., no filter applied).
        :return: A list of ARCSIScnProcess objects.

        """
        logger.debug("Creating Database Engine.")
        db_engine = sqlalchemy.create_engine(self.sqlite_db_conn, pool_pre_ping=True)
        logger.debug("Creating Database Session.")
        session_sqlalc = sqlalchemy.orm.sessionmaker(bind=db_engine)
        ses = session_sqlalc()
        logger.debug("Created Database Engine and Session.")

        logger.debug("Perform query to find scene.")
        query_result = (
            ses.query(ARCSIScnProcess)
            .filter(ARCSIScnProcess.geo_str_id == geo_str_id)
            .all()
        )
        scns = list()
        for scn in query_result:
            scns.append(scn)
        ses.close()
        logger.debug("Closed the database session.")
        return scns

    def set_scn_downloaded(self, product_id, sensor, download_path):
        """
        A function which sets that an the scene has been downloaded.

        :param product_id: the product ID for the scene
        :param sensor: the sensor for the scene
        :param download_path: the directory path for the downloaded data for the scene.

        """
        logger.debug("Creating Database Engine.")
        db_engine = sqlalchemy.create_engine(self.sqlite_db_conn, pool_pre_ping=True)
        logger.debug("Creating Database Session.")
        session_sqlalc = sqlalchemy.orm.sessionmaker(bind=db_engine)
        ses = session_sqlalc()
        logger.debug("Created Database Engine and Session.")

        logger.debug("Perform query to find scene.")
        query_result = (
            ses.query(ARCSIScnProcess)
            .filter(
                ARCSIScnProcess.product_id == product_id,
                ARCSIScnProcess.sensor == sensor,
            )
            .one_or_none()
        )
        if query_result is not None:
            query_result.download = True
            query_result.download_path = download_path
            ses.commit()
        ses.close()
        logger.debug("Closed the database session.")

    def get_scns_download(self, geo_str_id=None):
        """
        A function which gets a list of scenes which have not been downloaded.

        :param geo_str_id: Optionally, a geographic filter can be applied which is dependent on the sensor.
                           If landsat, this is 'r{row}_p{path}. If sentinel-2 then it is the granule.'
                           Default: None (i.e., no filter applied).
        :return: A list of ARCSIScnProcess objects.

        """
        logger.debug("Creating Database Engine.")
        db_engine = sqlalchemy.create_engine(self.sqlite_db_conn, pool_pre_ping=True)
        logger.debug("Creating Database Session.")
        session_sqlalc = sqlalchemy.orm.sessionmaker(bind=db_engine)
        ses = session_sqlalc()
        logger.debug("Created Database Engine and Session.")

        logger.debug("Perform query to find scene.")
        if geo_str_id is None:
            query_result = (
                ses.query(ARCSIScnProcess)
                .filter(ARCSIScnProcess.download == False)
                .all()
            )
        else:
            query_result = (
                ses.query(ARCSIScnProcess)
                .filter(
                    ARCSIScnProcess.geo_str_id == geo_str_id,
                    ARCSIScnProcess.download == False,
                )
                .all()
            )
        scns = list()
        for scn in query_result:
            scns.append(scn)
        ses.close()
        logger.debug("Closed the database session.")
        return scns

    def is_scn_downloaded(self, product_id, sensor):
        """
        A function to check whether a scene has been downloaded.

        :param product_id: the product ID for the scene
        :param sensor: the sensor for the scene
        :return: boolean (True: has been downloaded. False: has not been downloaded.

        """
        logger.debug("Creating Database Engine.")
        db_engine = sqlalchemy.create_engine(self.sqlite_db_conn, pool_pre_ping=True)
        logger.debug("Creating Database Session.")
        session_sqlalc = sqlalchemy.orm.sessionmaker(bind=db_engine)
        ses = session_sqlalc()
        logger.debug("Created Database Engine and Session.")

        logger.debug("Perform query to find scene.")
        query_result = (
            ses.query(ARCSIScnProcess)
            .filter(
                ARCSIScnProcess.product_id == product_id,
                ARCSIScnProcess.sensor == sensor,
            )
            .one_or_none()
        )
        ses.close()
        logger.debug("Closed the database session.")
        downloaded = False
        if query_result is not None:
            downloaded = query_result.download
        return downloaded

    def set_scn_ard(self, product_id, sensor, ard_path):
        """
        A function which sets that an ARD product has been generated for the scene.

        :param product_id: the product ID for the scene
        :param sensor: the sensor for the scene
        :param ard_path: the directory path for the ARD product.

        """
        logger.debug("Creating Database Engine.")
        db_engine = sqlalchemy.create_engine(self.sqlite_db_conn, pool_pre_ping=True)
        logger.debug("Creating Database Session.")
        session_sqlalc = sqlalchemy.orm.sessionmaker(bind=db_engine)
        ses = session_sqlalc()
        logger.debug("Created Database Engine and Session.")

        logger.debug("Perform query to find scene.")
        query_result = (
            ses.query(ARCSIScnProcess)
            .filter(
                ARCSIScnProcess.product_id == product_id,
                ARCSIScnProcess.sensor == sensor,
            )
            .one_or_none()
        )
        if query_result is not None:
            query_result.ard = True
            query_result.ard_path = ard_path
            ses.commit()
        ses.close()
        logger.debug("Closed the database session.")

    def get_scns_ard(self, geo_str_id=None):
        """
        A function which returns a list of scenes for which an ARD product needs to be generated.

        :param geo_str_id: Optionally, a geographic filter can be applied which is dependent on the sensor.
                           If landsat, this is 'r{row}_p{path}. If sentinel-2 then it is the granule.'
                           Default: None (i.e., no filter applied).
        :return: A list of ARCSIScnProcess objects.

        """
        logger.debug("Creating Database Engine.")
        db_engine = sqlalchemy.create_engine(self.sqlite_db_conn, pool_pre_ping=True)
        logger.debug("Creating Database Session.")
        session_sqlalc = sqlalchemy.orm.sessionmaker(bind=db_engine)
        ses = session_sqlalc()
        logger.debug("Created Database Engine and Session.")

        logger.debug("Perform query to find scene.")
        if geo_str_id is None:
            query_result = (
                ses.query(ARCSIScnProcess)
                .filter(ARCSIScnProcess.download == True, ARCSIScnProcess.ard == False)
                .all()
            )
        else:
            query_result = (
                ses.query(ARCSIScnProcess)
                .filter(
                    ARCSIScnProcess.geo_str_id == geo_str_id,
                    ARCSIScnProcess.download == True,
                    ARCSIScnProcess.ard == False,
                )
                .all()
            )
        scns = list()
        for scn in query_result:
            scns.append(scn)
        ses.close()
        logger.debug("Closed the database session.")
        return scns

    def get_processed_scns(self, geo_str_id=None):
        """
        A function which returns a list of scenes which have been processed (i.e., downloaded and ARD'd).

        :param geo_str_id: Optionally, a geographic filter can be applied which is dependent on the sensor.
                           If landsat, this is 'r{row}_p{path}. If sentinel-2 then it is the granule.'
                           Default: None (i.e., no filter applied).
        :return: A list of ARCSIScnProcess objects.

        """
        logger.debug("Creating Database Engine.")
        db_engine = sqlalchemy.create_engine(self.sqlite_db_conn, pool_pre_ping=True)
        logger.debug("Creating Database Session.")
        session_sqlalc = sqlalchemy.orm.sessionmaker(bind=db_engine)
        ses = session_sqlalc()
        logger.debug("Created Database Engine and Session.")

        logger.debug("Perform query to find scene.")
        if geo_str_id is None:
            query_result = (
                ses.query(ARCSIScnProcess)
                .filter(ARCSIScnProcess.download == True, ARCSIScnProcess.ard == True)
                .all()
            )
        else:
            query_result = (
                ses.query(ARCSIScnProcess)
                .filter(
                    ARCSIScnProcess.geo_str_id == geo_str_id,
                    ARCSIScnProcess.download == True,
                    ARCSIScnProcess.ard == True,
                )
                .all()
            )
        scns = list()
        for scn in query_result:
            scns.append(scn)
        ses.close()
        logger.debug("Closed the database session.")
        return scns

    def is_scn_ard(self, product_id, sensor):
        """
        A function to find whether an ARD product has been generated for the scene.

        :param product_id: the product ID for the scene
        :param sensor: the sensor for the scene
        :return: boolean (True: has been generated. False: has not been generated).

        """
        logger.debug("Creating Database Engine.")
        db_engine = sqlalchemy.create_engine(self.sqlite_db_conn, pool_pre_ping=True)
        logger.debug("Creating Database Session.")
        session_sqlalc = sqlalchemy.orm.sessionmaker(bind=db_engine)
        ses = session_sqlalc()
        logger.debug("Created Database Engine and Session.")

        logger.debug("Perform query to find scene.")
        query_result = (
            ses.query(ARCSIScnProcess)
            .filter(
                ARCSIScnProcess.product_id == product_id,
                ARCSIScnProcess.sensor == sensor,
            )
            .one_or_none()
        )
        ses.close()
        logger.debug("Closed the database session.")
        arded = False
        if query_result is not None:
            if query_result.ard_path is not None:
                arded = True
            else:
                arded = query_result.ard
        return arded

    def reset_all_scn(self, product_id, sensor, delpath=False):
        """
        A function to reset a scene (both ARD and download) - i.e., set it as not having been downloaded
        or converted to ARD.

        :param product_id: the product ID for the scene
        :param sensor: the sensor for the scene
        :param delpath: boolean; True: delete the download and ard paths if exists. Default: False

        """
        logger.debug("Creating Database Engine.")
        db_engine = sqlalchemy.create_engine(self.sqlite_db_conn, pool_pre_ping=True)
        logger.debug("Creating Database Session.")
        session_sqlalc = sqlalchemy.orm.sessionmaker(bind=db_engine)
        ses = session_sqlalc()
        logger.debug("Created Database Engine and Session.")

        logger.debug("Perform query to find scene.")
        query_result = (
            ses.query(ARCSIScnProcess)
            .filter(
                ARCSIScnProcess.product_id == product_id,
                ARCSIScnProcess.sensor == sensor,
            )
            .one_or_none()
        )
        if query_result is not None:
            logger.debug("Resetting ARD and Download fields for {}.".format(product_id))
            query_result.ard = False
            if (
                delpath
                and query_result.ard_path is not None
                and os.path.exists(query_result.ard_path)
            ):
                shutil.rmtree(query_result.ard_path)
            query_result.ard_path = None
            if (
                delpath
                and query_result.download_path is not None
                and os.path.exists(query_result.download_path)
            ):
                shutil.rmtree(query_result.download_path)
            query_result.download = False
            query_result.download_path = None
            ses.commit()
            logger.debug("Reset ARD and Download fields for {}.".format(product_id))
        else:
            logger.info("Failed to reset {} was not found.".format(product_id))
        ses.close()
        logger.debug("Closed the database session.")

    def reset_ard_scn(self, product_id, sensor, delpath=False):
        """
        A function to reset a scene as converted to ARD - i.e., set it as not having been converted to ARD.

        :param product_id: the product ID for the scene
        :param sensor: the sensor for the scene
        :param delpath: boolean; True: delete the ard path if exists. Default: False

        """
        logger.debug("Creating Database Engine.")
        db_engine = sqlalchemy.create_engine(self.sqlite_db_conn, pool_pre_ping=True)
        logger.debug("Creating Database Session.")
        session_sqlalc = sqlalchemy.orm.sessionmaker(bind=db_engine)
        ses = session_sqlalc()
        logger.debug("Created Database Engine and Session.")

        logger.debug("Perform query to find scene.")
        query_result = (
            ses.query(ARCSIScnProcess)
            .filter(
                ARCSIScnProcess.product_id == product_id,
                ARCSIScnProcess.sensor == sensor,
            )
            .one_or_none()
        )
        if query_result is not None:
            logger.debug("Resetting ARD fields for {}.".format(product_id))
            query_result.ard = False
            if (
                delpath
                and query_result.ard_path is not None
                and os.path.exists(query_result.ard_path)
            ):
                shutil.rmtree(query_result.ard_path)
            query_result.ard_path = None
            ses.commit()
            logger.debug("Reset ARD fields for {}.".format(product_id))
        else:
            logger.info("Failed to reset {} was not found.".format(product_id))
        ses.close()
        logger.debug("Closed the database session.")

    def reset_dwnld_scn(self, product_id, sensor, delpath=False):
        """
        A function which resets a download for a sensor - i.e., set it as not having been downloaded.

        :param product_id: the product ID for the scene
        :param sensor: the sensor for the scene
        :param delpath: boolean; True: delete the download path if exists. Default: False

        """
        logger.debug("Creating Database Engine.")
        db_engine = sqlalchemy.create_engine(self.sqlite_db_conn, pool_pre_ping=True)
        logger.debug("Creating Database Session.")
        session_sqlalc = sqlalchemy.orm.sessionmaker(bind=db_engine)
        ses = session_sqlalc()
        logger.debug("Created Database Engine and Session.")

        logger.debug("Perform query to find scene.")
        query_result = (
            ses.query(ARCSIScnProcess)
            .filter(
                ARCSIScnProcess.product_id == product_id,
                ARCSIScnProcess.sensor == sensor,
            )
            .one_or_none()
        )
        if query_result is not None:
            logger.debug("Resetting Download fields for {}.".format(product_id))
            query_result.download = False
            if (
                delpath
                and query_result.download_path is not None
                and os.path.exists(query_result.download_path)
            ):
                shutil.rmtree(query_result.download_path)
            query_result.download_path = None
            ses.commit()
            logger.debug("Reset Download fields for {}.".format(product_id))
        else:
            logger.info("Failed to reset {} was not found.".format(product_id))
        ses.close()
        logger.debug("Closed the database session.")
