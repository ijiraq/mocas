from unittest import TestCase

import numpy
import warnings
from mocas import core
from mp_ephem import EphemerisReader


class Test(TestCase):
    def test_search(self):
        warnings.simplefilter('ignore', category=UserWarning)
        ast_file = 'data/0123T02.ast'
        observations = EphemerisReader().read(ast_file)
        ssois_table = core.search(ast_file, start_date=observations[0].date, end_date=observations[-1].date)
        for observation in observations:
            frame = observation.comment.frame.split('p')[0]+'p'
            self.assertTrue(numpy.sum(ssois_table['Image'] == frame) > 0)
