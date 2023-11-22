import requests
from astropy.time import Time
from astropy.table import Table
from io import BytesIO
from mp_ephem import BKOrbit

class SSOIS:

    URL = "https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/cadcbin/ssos/ssosclf.pl"

    def __init__(self, obs_records, **kwargs):
        self.obs_records = obs_records
        self.orbit = BKOrbit(obs_records)
        self.params = {'obs': self.mpc_formatted_obs_record,
                       'search': 'bern',
                       'eunits': 'bern',
                       'extres': 'yes',
                       'xyres': 'yes',
                       'telinst': 'CFHT/MegaCam',
                       'format': 'tsv',
                       }
        self.params.update(kwargs)

    @property
    def mpc_formatted_obs_record(self):
        mpc_records = []
        for obs_record in self.obs_records:
            if obs_record.null_observation:
                continue
            mpc_records.append(obs_record.to_mpc()+"\r\n")
        return "".join(mpc_records)

    def __call__(self, start_date, end_date):
        self.params['epoch1'] = Time(start_date).isot
        self.params['epoch2'] = Time(end_date).isot
        response = requests.get(self.URL, params=self.params)
        response.raise_for_status()
        ssois_table = Table.read(BytesIO(response.content),
                                 format='ascii.csv',
                                 delimiter='\t',
                                 header_start=1,
                                 data_start=2)
        ssois_table['dra'] = 0.00
        ssois_table['ddec'] = 0.00
        ssois_table['PA'] = 0.00
        for row in ssois_table:
            self.orbit.predict(Time(row['MJD'], format='mjd'))
            row['dra'] = self.orbit.dra.value
            row['ddec'] = self.orbit.ddec.value
            row['PA'] = self.orbit.pa.value
        return ssois_table