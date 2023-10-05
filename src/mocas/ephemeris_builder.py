from astropy.coordinates import SkyCoord
import numpy
from astropy import units
from astropy.time import Time
from mp_ephem import BKOrbit

from .votable_file import TAPUploadVOTableFile


class SearchBoundsGenerator(object):
    """
    A class to produce an ephemeris  for a moving object.

    The format of the Ephemeris is list of two items [ [time_interval], [circle] ]

    where circle is centred on the mean coordinate of the sky path of the moving object and has a radius equal to the
    maximum separation of the sky path from the mean coordinate.

    """

    def __init__(self, orbit, start_date:Time, end_date:Time, time_step=1 * units.day, interval_size=30 * units.day):
        """
        Returns a time interval and an IVOA CIRCLE that is the union of sky locations visited by
        the SSO during that time interval.  The CIRCLE is centred on the mean coordinate of the sky path of the SSO
        and takes into account the uncertainty in the SSO's ephemeris position along a trajectory over time interval.

        Parameters
        ----------
        orbit: BKOrbit
            The orbit object of the SSO, used to compute the ephemeris.
        time_step : astropy.units.Quantity
            The minimum separation between two consecutive RA/Dec pairs in the ephemeris.
        """
        self.orbit = orbit
        self.time_step = time_step
        self.start_date = start_date
        self.end_date = end_date
        self.interval_size = interval_size
        self._coords = None
        self._interval_start_date = None
        self._interval_end_date = None

    @property
    def field_definitions(self) -> dict:
        return {'time_interval': {'datatype': 'double', 'arraysize': '2', 'xtype': 'interval'},
                'pos_circle': {'datatype': 'double', 'arraysize': '3', 'xtype': 'circle'}}

    @property
    def interval_start_date(self) -> Time:
        return self._interval_start_date

    @interval_start_date.setter
    def interval_start_date(self, value: Time):
        assert isinstance(value, Time)
        self._interval_start_date = value
        self._coords = None

    @property
    def interval_end_date(self) -> Time:
        return self._interval_end_date

    @interval_end_date.setter
    def interval_end_date(self, value: Time):
        assert isinstance(value, Time)
        self._interval_end_date = value
        self._coords = None

    @property
    def rows(self) -> list:
        """
        Return a list of interval rows for the ephemeris.
        :return:
        """
        rows = []
        self.interval_start_date = self.start_date
        while self.interval_start_date < self.end_date:
            self.interval_end_date = self.interval_start_date + self.interval_size
            rows.append(((self.interval_start_date.mjd,
                          self.interval_end_date.mjd),
                         self._circle_region))
            self.interval_start_date = self.interval_end_date
        return rows

    @property
    def _circle_region(self):
        """
        The centre of the circle is the mean coordinate of the sky path of the moving object.
        The radius of the circle is equal to the maximum separation of the sky path from the mean coordinate.

        Assumes that the set of sky coordinates does not have a mean location of 0,0,0 in cartesian coordinates.
        """
        coords = self._ephemeris_coordinates
        if len(coords) == 0:
            return None
        avg_coord = SkyCoord(coords.represent_as('cartesian').mean(),
                             representation_type='spherical')
        radius = avg_coord.separation(coords).max().to('degree').value
        return avg_coord.ra.degree, avg_coord.dec.degree, radius

    @property
    def _ephemeris_coordinates(self) -> SkyCoord:
        """
        Return a SkyCoord list of the RA/Dec pairs of the SSO locations during the time interval.

        Treats the ra/dec uncertainty as a box and computes the corners of the box for each time step.
        """

        if self._coords is not None:
            return self._coords

        _coords = []
        box_corners = numpy.array([[-1, -1], [-1, 1], [1, 1], [1, -1]])
        for current_time in numpy.arange(self.interval_start_date.jd, self.interval_end_date.jd, self.time_step.to(units.day).value):
            self.orbit.predict(current_time)
            uncertainty_box_offsets = (box_corners *
                                       (self.orbit.dra.value, self.orbit.ddec.value) *
                                       (self.orbit.dra.unit, self.orbit.ddec.unit)).T
            _coords.extend(self.orbit.coordinate.spherical_offsets_by(uncertainty_box_offsets[0],
                                                                      uncertainty_box_offsets[1]))
        self._coords = SkyCoord(_coords)
        return self._coords

    @property
    def votable(self):
        """
        Return a VOTable representation of the ephemeris.
        """
        return TAPUploadVOTableFile(self.rows, self.field_definitions)
