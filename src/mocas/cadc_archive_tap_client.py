from io import BytesIO

from astropy.table import Table
from cadctap import CadcTapClient
import pyvo
import cadcutils


class CADCArchiveTapClient(CadcTapClient):
    """
    Convenience class to create a CadcTapClient for the CADC Archive TAP service accessed via a CERT
    """

    def __init__(self, subject=None):
        if subject is None:
            subject = cadcutils.net.Subject()
        super().__init__(resource_id=self.cadc_archive_resource_id,
                         subject=subject)

    @property
    def cadc_archive_resource_id(self) -> str:
        return pyvo.registry.search(keywords=["CADC"],
                                    servicetype="tap")["CADC TAP"].ivoid

    def get_table(self, query: str) -> pyvo.dal.tap.TAPResults:
        """
        Get a table from the CADC Archive TAP service.
        """
        with BytesIO() as f:
            self.query(query, output_file=f)
            f.seek(0)
            table = Table.read(f, format='votable')
        return table
