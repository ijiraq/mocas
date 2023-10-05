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
