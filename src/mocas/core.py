"""
Do a search of the CADC CAOM database for a moving object.
"""
import argparse
import glob
import logging
import os
import time
import warnings
from io import BytesIO
from tempfile import NamedTemporaryFile

from astropy import wcs
from astropy.table import Table, unique
from astropy.time import Time
from mp_ephem import BKOrbit

from .cadc_archive_tap_client import CADCArchiveTapClient
from .ephemeris_builder import SearchBoundsGenerator
from .votable_file import TAPUploadVOTableFile


def get_orbit_from_ast_file(ast_file) -> BKOrbit:  # noqa: N802
    """
    Fit an orbit to a set of observations in an ast_file and return as object.
    :param ast_file:
    :return:
    """
    return BKOrbit(None, ast_filename=ast_file)


def query_cadc_archive_against_votable(query: str, votable_file: TAPUploadVOTableFile) -> Table:
    """
    Run QUERY against CADC ARCHIVE TAP service using the TAPUploadVOTableFile as the tmptable.
    :param query:
    :param votable_file:
    :return:
    """
    with NamedTemporaryFile() as f:
        votable_file.to_xml(f)
        f.flush()
        result = BytesIO()
        CADCArchiveTapClient().query(query=query,
                                     output_file=result,
                                     tmptable=f"tmptable:{f.name}",
                                     timeout=10)
    result.seek(0)
    return Table.read(result, format='votable')


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
        "SELECT observationID, target_name, (time_bounds_lower+time_bounds_upper)/2 as mid_mjd "
        "FROM caom2.Plane AS p JOIN caom2.Observation AS o ON p.obsID=o.obsID "
        ", tap_upload.tmptable AS arc "
        "WHERE o.instrument_name = '{}' "
        "AND o.collection = '{}' "
        "AND INTERSECTS( arc.time_interval, p.time_bounds_samples ) = 1  "
        "AND INTERSECTS( arc.pos_circle, p.position_bounds) = 1")
    result = query_cadc_archive_against_votable(
        _search_along_ephemeris_query.format(instrument, collection),
        votable_file)
    if len(result) > 0:
        result = unique(result, keys='observationID')
    return result


def get_observation_details(observation_ids: list) -> Table:
    """
    Get the WCS of a set of observations
    :param observation_ids:
    :return:
    """
    _observation_details_query = ("SELECT "
                                  "o.observationID, "
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
                                  "FROM "
                                  "caom2.Observation AS o "
                                  "JOIN caom2.Plane as p ON p.obsID=o.obsID "
                                  "JOIN caom2.Artifact AS a ON p.planeID=a.planeID "
                                  "JOIN caom2.Part AS pa ON a.artifactID=pa.artifactID "
                                  "JOIN caom2.Chunk AS c ON pa.partID=c.partID "
                                  ", tap_upload.tmptable AS obs "
                                  "WHERE o.observationID = obs.observationID AND p.calibrationLevel = 2 "
                                  "AND a.productType = 'science' ")
    field_definitions = {'observationID': {'datatype': 'char', 'arraysize': '*'}}
    observation_id_votable_file = TAPUploadVOTableFile(observation_ids, field_definitions)
    return query_cadc_archive_against_votable(_observation_details_query, observation_id_votable_file)


def filter_for_artifacts_containing_arc(artifacts: Table, orbit: BKOrbit) -> Table:
    """
    Filter through the list of observations (using the wcs of each artifact) to determine which contain the arc.

    :param artifacts:
    :param orbit:
    :return:
    """
    result = []
    for artifact in artifacts:
        orbit.predict(Time(artifact['mid_mjd'], format='mjd'))
        result.append(wcs.WCS(artifact).footprint_contains(orbit.coordinate))
    return artifacts[result]


ssois_column_names = ['Image', 'Filter', 'Image_Target',
                      'Ext', 'X', 'Y', 'MJD',
                      'Object_RA', 'Object_Dec', 'dra', 'ddec']


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
        ssois_columns['Image'].append(f"{artifact['observationID']}p")
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
    return Table(ssois_columns)


def search(ast_file: str, start_date: Time, end_date: Time):
    """
    Search the CADC archive for observations of a moving object in a time interval.

    :param ast_file:
    :param start_date:
    :param end_date:
    :return:
    """
    orbit = get_orbit_from_ast_file(ast_file)
    observations = search_along_arc(SearchBoundsGenerator(orbit, start_date, end_date).votable)
    artifact_wcs = get_observation_details(list(observations["observationID"]))
    return build_ssois_table(filter_for_artifacts_containing_arc(artifact_wcs, orbit),
                             orbit)


def main():
    parser = argparse.ArgumentParser(description='Search the CADC CFHT MegaPrime archive for observations of '
                                                 'a moving object in a time interval.')
    parser.add_argument('ast_file', type=str, help='The path to the ast file containing '
                                                   'the observations of the ')
    parser.add_argument('start_date', type=str, help='The start date of the time interval to search')
    parser.add_argument('end_date', type=str, help='The end date of the time interval to search')
    parser.add_argument('--output', type=str, help='The path to the output file')
    parser.add_argument('--log-level', type=str, help='The log level to use',
                        default='INFO',
                        choices=['DEBUG',
                                 'INFO',
                                 'WARNING',
                                 'ERROR',
                                 'CRITICAL'])
    args = parser.parse_args()
    warnings.simplefilter('ignore', category=UserWarning)
    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s %(levelname)s %(message)s')
    start_date = Time(args.start_date)
    end_date = Time(args.end_date)
    ast_file = args.ast_file
    start_time = time.time()
    logging.info(f"Searching for observations along path {ast_file}")
    search(ast_file, start_date, end_date).write(f"{ast_file.rstrip('ast')}tsv",
                                                 format='ascii',
                                                 delimiter='\t',
                                                 overwrite=True)
    end_time = time.time()
    logging.info(f"Time taken: {end_time - start_time}")


if __name__ == "__main__":
    main()
