import numpy as np
import astropy.units as u
from scipy import special

class Event(object):

    def __init__(
    self,
    fluence,
    f0,
    wf,
    freqs
    ):


        self.f0 = f0
        self.wf = wf
        self.fluence = fluence
        self.freqs = freqs

    def get_spec(self):

        """
        Get FRB spectrum as a toy Gaussian function
        """

        f  = self.freqs
        Wf = self.wf
        f0 = self.f0
        F  = self.fluence

        sigma = Wf.to(u.MHz).value / 2 / np.sqrt(2 * np.log(2))
        f = f.to(u.MHz).value
        f0 = f0.to(u.MHz).value
        F = F.to(u.Jy * u.ms).value

        spec = (F / (np.sqrt(2 * np.pi) * sigma)) * np.exp(
            -(1 / 2) * ((f - f0) / sigma) ** 2
        )

        return spec

    def fluence_perband(self, fstart, fend):

        """
        Get the fluence contribution in a given sub-band
        """

        fstart = fstart.to(u.MHz).value
        fend   = fend.to(u.MHz).value

        mu_f = self.f0
        Wf = self.wf
        F = self.fluence


        mu_f = mu_f.to(u.MHz).value
        sig_f = Wf.to(u.MHz).value / 2 / np.sqrt(2 * np.log(2))
        F = F.value

        inp_start = (fstart - mu_f) / (np.sqrt(2 * np.pi) * sig_f)
        inp_end = (fend - mu_f) / (np.sqrt(2 * np.pi) * sig_f)
        integral = 0.5 * (special.erf(inp_end) - special.erf(inp_start))

        Fmeasured = F *  integral

        return Fmeasured

    def plot_event(self):

        """
        Make a simple plot 
        """

        freqs = self.freqs

        spec = self.get_spec()

        plt.figure(figsize = (8,4))
        plt.ylabel("Fluence (a.u.)")
        plt.xlabel("Frequency (MHz)")
        plt.plot(freqs.value, spec)
        plt.show()
