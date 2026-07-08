"""
This project originally featured a section where we analyzed our entire dataset
and grouped into galactic components for the sake of a visualization of where 
different groups sat on a CMD. However, this proceeded in a very 
disorganized fashion with multiple rounds (as it was an earlier project in my coding 
journey), which is inefficient. So we utilized object-oriented programming to solve this. 

Additionally, assessing the performance after each selection was important in the initial
project, so the "rounds" will be kept to allow this method to remain as a feature.

Basically, a "reusable tool" for selecting cluster members from Gaia
archive using sky position and proper-motion cuts.

Additionally, can determine which galactic component the cluster lives in.

(Update: this now includes optional existing radvels and well-measured parallax cuts!)
(Update 2: included optional additional masking, like for red clump! In order to use this, it must be a boolean array of the same length.)

"""

import numpy as np
import astropy.coordinates as coord
import astropy.units as u
from astropy.table import Table, QTable

class SelectionRound:
	"""
	This class is responsible for creating/holding new subsets of stars after a selection round.
	New instances become new subsets, based on parameters of select_round method within ClusterSelector class.
	Also can analyze the subset, but cannot modify it.

	"""

	def __init__(self, round_num, ra, dec, mag, color, parallax, parallax_er, pmra, pmdec, pmra_er, pmdec_er, radvel, radvel_er, teff):
		self.round_num = round_num
		self.ra = ra
		self.dec = dec
		self.mag = mag
		self.color = color
		self.parallax = parallax
		self.parallax_er = parallax_er
		self.pmra = pmra
		self.pmdec = pmdec
		self.pmra_er = pmra_er
		self.pmdec_er = pmdec_er
		self.radvel = radvel
		self.radvel_er = radvel_er
		self.teff = teff

	def __len__(self): # included for debugging purposes
		return len(self.ra)
		

	def galactic_comp(self):
		dist = coord.Distance(parallax = u.Quantity(self.parallax))
		c = coord.SkyCoord(ra = self.ra, dec = self.dec, distance = dist, pm_ra_cosdec = self.pmra, pm_dec = self.pmdec, radial_velocity = self.radvel)
		
		coord.Galactocentric()
		
		galcen = c.transform_to(coord.Galactocentric(z_sun=15.3*u.pc,
											galcen_distance=8.18*u.kpc,
											galcen_v_sun=coord.CartesianDifferential(np.array([11.1,250.1,7.25]),unit=u.kilometer/u.second)))
		R = np.sqrt(galcen.x.value**2 + galcen.y.value**2)				 

		thin_disk = (np.abs(galcen.z.value) < 200)
		thick_disk = (np.abs(galcen.z.value) >= 200) & (np.abs(galcen.z.value) < 1000)
		halo = (np.abs(galcen.z.value) >= 1000) & (np.abs(galcen.z.value) < 30000)

		labels = np.select([thin_disk, thick_disk, halo],["thin_disk", "thick_disk", "halo"],default="unclassified")
		return labels, dist, galcen, R
		

class ClusterSelector:
	"""
	This holds the full dataset initially, then gets passed to the SelectionRound class.
	Is the tool that will perform the selections.

	Each call to the function 'select_round' applies a cut based on sky position and
	(optionally) a proper-motion window centered on the mean pmra/pmdec of
	a previous round (sigma-clipped). Reduces repetitive code of previous version.

	"""

	def __init__(self, data):
		if isinstance(data, (dict, Table, QTable)):
			self.rounds = []
			self.current_round = 0
			self.ra = data["ra"]
			self.dec = data["dec"]
			self.mag = data["phot_g_mean_mag"]
			self.color = data["bp_rp"]
			self.parallax = data["parallax"]
			self.parallax_er = data["parallax_error"]
			self.pmra = data["pmra"]
			self.pmdec = data["pmdec"]
			self.pmra_er = data["pmra_error"]
			self.pmdec_er = data["pmdec_error"]
			self.radvel = data["radial_velocity"]
			self.radvel_er = data["radial_velocity_error"]
			self.teff = data["teff_gspphot"] if "teff_gspphot" in data.colnames else None

		elif isinstance(data, SelectionRound):
			self.rounds = [data]
			self.current_round = data.round_num
			self.ra = data.ra
			self.dec = data.dec
			self.mag = data.mag
			self.color = data.color
			self.parallax = data.parallax
			self.parallax_er = data.parallax_er
			self.pmra = data.pmra
			self.pmdec = data.pmdec
			self.pmra_er = data.pmra_er
			self.pmdec_er = data.pmdec_er
			self.radvel = data.radvel
			self.radvel_er = data.radvel_er
			self.teff = data.teff
		else:
			raise TypeError("ClusterSelector expects Astropy Table, QTable, dict, or SelectionRound data types. Modify class or use different type please!")

	# sigma = std dev
	# ra_range and dec_range should be arrays

	def select(self, ra_range=(0, 360), dec_range=(-90, 90), ref_round=None, sigma_window=None, radvel_mask=False, par_mask=False, add_mask=None):
		indices = (self.ra >= ra_range[0]) & (self.ra <= ra_range[1]) & (self.dec >= dec_range[0]) & (self.dec <= dec_range[1])

		# determine if proper motion clipping will be used
		if sigma_window is not None:
			if ref_round is not None:
				if sigma_window[0] == 'pm':
					mean_pmra = np.mean(ref_round.pmra)
					mean_pmdec = np.mean(ref_round.pmdec)

					std_pmra = np.std(ref_round.pmra) * sigma_window[1]
					std_pmdec = np.std(ref_round.pmdec) * sigma_window[1]

					upper_pmra = mean_pmra + std_pmra
					lower_pmra = mean_pmra - std_pmra
					upper_pmdec = mean_pmdec + std_pmdec
					lower_pmdec = mean_pmdec - std_pmdec
					
					indices &=  ((self.pmra >= lower_pmra) & (self.pmra <= upper_pmra)) & ((self.pmdec >= lower_pmdec) & (self.pmdec <= upper_pmdec))
				elif sigma_window[0] == 'color':
					print("In progress")
			else:
				raise ValueError("Where is your previous round's data? Don't run sigma clipping on first round")

		# typically will be used to determine where radial velocity exists
		if radvel_mask:
			indices &= ~self.radvel.mask

		# to choose well-measured parallax
		if par_mask:
			valid_parallax = (np.isfinite(self.parallax) & np.isfinite(self.parallax_er) & (self.parallax_er > 0) & (self.parallax > 0))

			indices &= valid_parallax
			indices &= ((self.parallax / self.parallax_er) > 10)
		if add_mask is not None:
			indices &= add_mask

		idx = np.where(indices)

		result = SelectionRound(round_num = self.current_round + 1,ra = self.ra[idx], dec = self.dec[idx], mag = self.mag[idx], color = self.color[idx], parallax = self.parallax[idx],
			parallax_er = self.parallax_er[idx], pmra = self.pmra[idx], pmdec = self.pmdec[idx], pmra_er = self.pmra_er[idx], pmdec_er = self.pmdec_er[idx],
			radvel = self.radvel[idx], radvel_er = self.radvel_er[idx], teff = self.teff[idx] if self.teff is not None else None)
		
		self.current_round = result.round_num
		self.rounds.append(result)

		return result


def absolute_magnitude(parallax, apparent_mag):

	mean_parallax = np.mean(parallax)
	mean_distance_pc = 1000 / mean_parallax
	return apparent_mag - 5 * (np.log10(mean_distance_pc) - 1)
