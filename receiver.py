import numpy as np
import astropy.units as u

class Receiver(object):

    def __init__(
    self,
    name,
    fc,
    bw,
    nchan,
    dt,
    gain,
    tsys,
    npol
    ):

        self.name = name
        self.fc = fc
        self.bw = bw
        self.nchan = nchan
        self.tsys = tsys
        self.gain = gain
        self.dt = dt
        self.npol = npol

    def get_freqs(self):

        """
        Get frequency array
        """

        freqs = np.linspace(self.fc - self.bw / 2 , self.fc + self.bw / 2, self.nchan)

        return freqs


    def get_mdf(self, threshold = 10, width = 1 * u.ms):

        """
        Get Mimimum Detectable Fluence in Jy ms
        """

        Fmin  = threshold * self.tsys / np.sqrt(self.npol * self.bw * width) / self.gain * width

        return Fmin.to(u.Jy * u.ms)
