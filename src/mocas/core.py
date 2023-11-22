"""
Do a search of the CADC CAOM database for a moving object.
"""
import asyncio
import logging
import os
import time
from io import BytesIO
from tempfile import NamedTemporaryFile
from astropy import wcs
from astropy.table import Table, unique, vstack
from astropy.time import Time
from mp_ephem import BKOrbit, EphemerisReader
from .cadc_archive_tap_client import CADCArchiveTapClient
from .ephemeris_builder import SearchBoundsGenerator
from .votable_file import TAPUploadVOTableFile
from .ssois import SSOIS

def get_orbit_from_ast_file(ast_file) -> BKOrbit:  # noqa: N802
    """
    Fit an orbit to a set of observations in an ast_file and return as object.
    :param ast_file:
    :return:
    """
    return BKOrbit(None, ast_filename=ast_file)


async def query_cadc_archive_against_votable(query: str,
                                             votable_file: TAPUploadVOTableFile,
                                             query_name='query',
                                             chunk_size=20) -> Table:
    """
    Run QUERY against CADC ARCHIVE TAP service using the TAPUploadVOTableFile as the tmptable.
    :param query:
    :param votable_file:
    :param query_name: string to use in logging
    :param chunk_size: number of rows to send in each query
    :return:
    """
    arc_length = len(votable_file)
    logging.info(f"{query_name} query has {arc_length} rows")
    chunks = arc_length//chunk_size + 1
    logging.info(f"Breaking into {chunks} asynchronous queries")
    vo_tables = []
    for chunk in range(chunks):
        vo_tables.append(TAPUploadVOTableFile(votable_file.get_first_table().array[chunk:chunk + chunk_size],
                                              votable_file.field_definitions))
    logging.info(f"Got {len(vo_tables)} sets")
    return vstack(await asyncio.gather(*[async_query_cadc_archive_against_votable(query, vo_table)
                                         for vo_table in vo_tables]))


async def async_query_cadc_archive_against_votable(query: str, votable_file: TAPUploadVOTableFile) -> Table:
    """
    """
    client = CADCArchiveTapClient()
    logging.info(f"Querying against votable of length: {len(votable_file)}")
    with NamedTemporaryFile(delete=False) as f:
        votable_file.to_xml(f)
        start_of_query = time.time()

    with NamedTemporaryFile() as result:
        logging.info(f"Starting query")
        await asyncio.to_thread(client.query,
                                query=query,
                                output_file=result.name,
                                tmptable=f"tmptable:{f.name}",
                                timeout=60)
        end_of_query = time.time()
        logging.info(f"Query took {end_of_query - start_of_query} seconds")
        result.seek(0)
        result_table = Table.read(result.name, format='votable')
    os.remove(f.name)
    logging.debug(f"{result_table}")
    return result_table


def search_along_arc(votable_file: TAPUploadVOTableFile,
                     collection: str = "CFHT",
                     instrument: str = "MegaPrime") -> Table:
    """
    Search the CADC archive for observations of a moving object in a time interval.

    :param votable_file: The arc to search along in VOTable format
    :param instrument: The instrument to search within
    :param collection: The collection (telescope data) to search within
    :return:
    """
    _search_along_ephemeris_query = (
        "SELECT "
        " planeID "
        "FROM caom2.Plane AS p "
        ", tap_upload.tmptable AS arc "
        "WHERE "
        "p.calibrationLevel = 2 "
        "AND INTERSECTS( arc.time_interval, p.time_bounds_samples ) = 1  "
        "AND INTERSECTS( arc.polygon, p.position_bounds) = 1")
    result = asyncio.run(query_cadc_archive_against_votable(
        _search_along_ephemeris_query.format(instrument, collection),
        votable_file, query_name='search_along_ephemeris'))
    if len(result) > 0:
        result = unique(result, keys='planeID')
    logging.info(f"Found {len(result)} observation planes overlapping the ephemeris")
    return result


def get_artifact_wcs(plane_ids: list) -> Table:
    """
    Get the WCS of a set of observations
    :param plane_ids:
    :return:
    """
    _observation_details_query = ("SELECT "
                                  "o.observationID as observation_id,"
                                  "p.productID as Image, "
                                  "o.target_name as Image_Target, "
                                  "p.energy_bandpassName as Filter, "
                                  "(p.time_bounds_upper + p.time_bounds_lower)/2 as mid_mjd, "
                                  "pa.name as EXTNAME, "
                                  "position_coordsys AS COORDSYS, "
                                  "position_equinox AS EQUINOX, "
                                  "position_axis_function_dimension_naxis1 AS NAXIS1, "
                                  "position_axis_axis1_ctype AS CTYPE1, "
                                  "position_axis_axis1_cunit AS CUNIT1, "
                                  "position_axis_function_refCoord_coord1_pix AS CRPIX1, "
                                  "position_axis_function_refCoord_coord1_val AS CRVAL1, "
                                  "position_axis_function_dimension_naxis2 AS NAXIS2, "
                                  "position_axis_axis2_ctype AS CTYPE2, "
                                  "position_axis_axis2_cunit AS CUNIT2, "
                                  "position_axis_function_refCoord_coord2_pix AS CRPIX2, "
                                  "position_axis_function_refCoord_coord2_val AS CRVAL2, "
                                  "position_axis_function_cd11 as CD1_1, "
                                  "position_axis_function_cd12 as CD1_2, "
                                  "position_axis_function_cd21 as CD2_1, "
                                  "position_axis_function_cd22 as CD2_2 "
                                  "FROM caom2.Plane AS p "
                                  "JOIN caom2.Observation AS o ON o.obsID=p.obsID "
                                  "JOIN caom2.Artifact AS a ON p.planeID=a.planeID "
                                  "JOIN caom2.Part AS pa ON a.artifactID=pa.artifactID "
                                  "JOIN caom2.Chunk AS c ON pa.partID=c.partID "
                                  ", tap_upload.tmptable AS obs "
                                  "WHERE p.planeID = obs.eph_planeID "
                                  "AND a.productType = 'science' ")
    # obsID_list = "("+",".join(["'"+plane_id+"'" for plane_id in plane_ids])+")"
    # _observation_details_query = _observation_details_query.format(obsID_list)
    field_definitions = {'eph_planeID': {'datatype': 'char', 'arraysize': '36', 'xtype': 'uuid'}}
    observation_id_votable_file = TAPUploadVOTableFile(plane_ids, field_definitions)
    # return query_cadc_archive(_observation_details_query)
    return asyncio.run(query_cadc_archive_against_votable(_observation_details_query, observation_id_votable_file,
                                                          query_name='observation_details'))


def filter_for_artifacts_containing_arc(artifacts_wcs: Table, orbit: BKOrbit, dblist: list) -> Table:
    """
    Filter through the list of observations (using the wcs of each artifact) to determine which contain the arc.

    :param artifacts_wcs:
    :param orbit:
    :param dblist:  a list of observation_id values in the dbimages directory.
    :return:
    """
    result = []
    for artifact_wcs in artifacts_wcs:
        result.append(False)
        if dblist is not None and artifact_wcs['observation_id'] not in dblist:
            logging.warning(f"Skipping {artifact_wcs['observation_id']} as not in observation id list.")
            continue
        try:
            orbit.predict(Time(artifact_wcs['mid_mjd'], format='mjd'))
            this_artifact_wcs = wcs.WCS(artifact_wcs)
            x, y = this_artifact_wcs.wcs_world2pix(orbit.coordinate.ra.degree, orbit.coordinate.dec.degree, 1)
            if 32 < x < 2080 and 0 < y < 4612:
                result[-1] = True
        except Exception as e:
            logging.warning(f"Failed to build a wcs using {artifact_wcs}: {e}")
    if sum(result) > 0:
        return unique(artifacts_wcs[result], keys=["Image", "Filter", "Image_Target", "EXTNAME"])
    else:
        return artifacts_wcs[result]


ssois_column_names = ['Image', 'Filter', 'Image_Target',
                      'Ext', 'X', 'Y', 'MJD',
                      'Object_RA', 'Object_Dec', 'dra', 'ddec', 'PA']


def build_ssois_table(artifacts: Table, orbit: BKOrbit) -> Table:
    """
    construct the SSOIS table from the orbit and the observation and wcs info for each artifact.
    :param orbit:
    :param artifacts:
    :return:
    """
    ssois_columns = dict([(name, []) for name in ssois_column_names])
    for artifact in artifacts:
        orbit.predict(Time(artifact['mid_mjd'], format='mjd'))
        x, y = wcs.WCS(artifact).wcs_world2pix(orbit.coordinate.ra.degree, orbit.coordinate.dec.degree, 1)
        ssois_columns['Image'].append(f"{artifact['Image']}")
        ssois_columns['Filter'].append(artifact['Filter'])
        ssois_columns['Image_Target'].append(artifact['Image_Target'])
        ssois_columns['Ext'].append(artifact['EXTNAME'])
        ssois_columns['X'].append(x)
        ssois_columns['Y'].append(y)
        ssois_columns['MJD'].append(artifact['mid_mjd'])
        ssois_columns['Object_RA'].append(orbit.coordinate.ra.degree)
        ssois_columns['Object_Dec'].append(orbit.coordinate.dec.degree)
        ssois_columns['dra'].append(orbit.dra.value)
        ssois_columns['ddec'].append(orbit.ddec.value)
        ssois_columns['PA'].append(orbit.pa.value)
    return Table(ssois_columns)


def search(ast_file: str, start_date: Time, end_date: Time, observation_ids: list = None, use_ssois=False) -> Table:
    """
    Search the CADC archive for observations of a moving object in a time interval.

    :param ast_file:
    :param start_date:
    :param end_date:
    :return:
    """
    if use_ssois:
        ssois_search = SSOIS(EphemerisReader().read(ast_file))
        return ssois_search(start_date, end_date)
    orbit = get_orbit_from_ast_file(ast_file)
    observations = search_along_arc(SearchBoundsGenerator(orbit, start_date, end_date).votable)
    logging.info(f"Found {len(observations)} possible observations, checking if on detector.")
    artifact_wcs = get_artifact_wcs(list(observations["planeID"]))
    logging.info(f"Found {len(artifact_wcs)} extensions to check.")
    return build_ssois_table(filter_for_artifacts_containing_arc(artifact_wcs, orbit, observation_ids),
                             orbit)
