import numpy as np 
import matplotlib
import ternary 
from scipy.interpolate import interp1d
from powerreader import * 
from numpy.random import rand 

#Set up plot style
font = {'size'   : 20}
matplotlib.rc('font', **font)
matplotlib.rc('font', serif='Computer Modern Roman')

#Define colours 
martaRed = "#c24c51" 
martaGreen = "#54a666"
martaBlue = "#4c70b0"
martaPurple = "#7f70b0"
martaGold = "#ccb873"
Icolour = '#DB2420' #Central line red
Ocolour = '#00A0E2' #Victoria line blue
O2colour = '#868F98' #Jubilee line grey
D2colour = '#F386A0' #Hammersmith line pink 
D4colour = '#97015E' #Metropolitan line magenta
D6colour = '#B05F0F' #Bakerloo line brown
D5colour = '#00843D' #District line green


class Battery: 

	def __init__(self, R, C, **kwargs): 
		self.__dict__.update(kwargs) 
		self.max_rate = R 
		self.max_capacity = C 
		self.charge_level = 0.0 

	def charge_battery(self, p, **kwargs):
		self.__dict__.update(kwargs)
		#Check if this is more than max rate 
		charge_before = self.charge_level
		if p > self.max_rate: 
			#Then charge at the max rate, if there is space
			x = self.charge_level + self.max_rate
			if x >= self.max_capacity: 
				self.charge_level = self.max_capacity 
			else:
				self.charge_level = x
		else:
			#Then take everything in
			x = self.charge_level + p
			if x >= self.max_capacity: 
				self.charge_level = self.max_capacity 
			else:
				self.charge_level = x
		self.delta_charge = self.charge_level - charge_before 
		return p - self.delta_charge 

	def discharge_battery(self, p, **kwargs):
		self.__dict__.update(kwargs) 
		charge_before = self.charge_level 
		if abs(p) > self.max_rate:
			x = self.charge_level - self.max_rate 
			if x <= 0.0:
				self.charge_level = 0.0 
			else:
				self.charge_level = x 
		else:
			x = self.charge_level - abs(p) 
			if x <= 0.0: 
				self.charge_level = 0.0
			else:
				self.charge_level = x 
		self.delta_charge = self.charge_level - charge_before 
		return p - self.delta_charge 


class House:

	def __init__(self, **kwargs):
		self.__dict__.update(kwargs)
		self.powertimes = [] 
		self.powervals = [] 
		self.generationtimes = [] 
		self.generationvals = [] 
		self.has_battery = False

	def make_power_interpolator(self, **kwargs):
		self.__dict__.update(kwargs)
		self.power_interpolator = interp1d(self.powertimes, self.powervals)

	def make_generation_interpolator(self, **kwargs):
		self.__dict__.update(kwargs) 
		self.generation_interpolator = interp1d(self.generationtimes, self.generationvals)

	def make_battery(self, R, C, **kwargs):
		self.__dict__.update(kwargs) 
		self.has_battery = True 
		self.powerwall = Battery(R, C) 

	def delete_battery(self, **kwargs):
		self.__dict__.update(kwargs)
		if self.has_battery:
			del self.powerwall
			self.has_battery = False

	def fill_missing_generation(self, t, **kwargs):
		self.__dict__.update(kwargs) 
		indeces_before_t = np.argwhere(self.generationtimes < t)
		try:
			generation_before_t = self.generationvals[indeces_before_t[-1]] 
		except IndexError: 
			generation_before_t = self.generationvals[0]
		return generation_before_t[0] 

	def get_states(self, t, **kwargs):
		self.__dict__.update(kwargs) 
		demand = self.power_interpolator(t) 
		try: 
			production = self.generation_interpolator(t)
		except ValueError:
			production = self.fill_missing_generation(t) 
		except AttributeError:
			production = 0.0
		storage = self.powerwall.delta_charge
		return demand, production, storage

	def get_house_power(self, t, **kwargs):
		self.__dict__.update(kwargs) 
		try: 
			x = self.generation_interpolator(t) - self.power_interpolator(t)
		except AttributeError:
			x = -self.power_interpolator(t)
		except ValueError:
			x = self.fill_missing_generation(t) - self.power_interpolator(t)
		if self.has_battery:
			if x > 0.0:
				x = self.powerwall.charge_battery(x)
			if x < 0.0: 
				x = self.powerwall.discharge_battery(x)
			return x
		else:
			return x 



class MicroGrid:   

	def __init__(self, **kwargs):
		self.__dict__.update(kwargs)
		self.n = 50 
		self.month = 1
		self.penetration = 49
		self.houses = []
		self.batteries_on_grid = False 

	def make_houses(self, **kwargs):
		self.__dict__.update(kwargs) 
		self.Pvec = np.zeros(self.n) 
		for i in range(self.n-1):
			h = House()
			self.houses.append(h) 

	def assign_powerprofiles(self, **kwargs):
		self.__dict__.update(kwargs)
		print("Assigning power profiles ... ")
		i = 0  
		for h in self.houses: 
			print(i)
			assigned = False 
			while assigned == False: 
				df, flag = make_random_week_profiles(self.month)
				assigned = flag 
			h.powertimes = np.array(df.secs) 
			h.powervals = np.array(df.means) 
			h.make_power_interpolator() 
			i+=1

	def assign_genprofiles(self, **kwargs):
		self.__dict__.update(kwargs) 
		i = 0
		for h in self.houses[:self.penetration]:
			print(i)
			assigned = False
			while assigned == False: 
				df, flag = make_random_week_profiles_PV(self.month)
				assigned = flag 
			h.generationtimes = np.array(df.secs) 
			h.generationvals = np.array(df.means)  
			h.make_generation_interpolator()
			i += 1

	def assign_batteries(self, **kwargs):
		self.__dict__.update(kwargs)
		self.batteries_on_grid = True 
		for h in self.houses[:self.penetration]:
			h.make_battery(2.5, 13.5)

	def get_av_battery_charge(self, **kwargs):
		self.__dict__.update(kwargs)
		if self.batteries_on_grid:
			charges = [] 
			for h in self.houses[:self.penetration]:
				bat = h.powerwall
				c = bat.charge_level 
				charges.append(c)
			return np.mean(charges) 
		else:
			return 0.0 

	def get_power_breakdown(self, t, **kwargs): 
		self.__dict__.update(kwargs) 
		self.gen_tot = 0 
		self.pow_tot = 0 
		for h in self.houses: 
			self.pow_tot += h.power_interpolator(t) 
		for h in self.houses[:self.penetration]: 
			try: 
				self.gen_tot += h.generation_interpolator(t)
			except ValueError:
				self.gen_tot += h.fill_missing_generation(t) 

	def get_hists(self, t, **kwargs):
		self.__dict__.update(kwargs) 
		self.gen_stats = [] 
		self.pow_stats = [] 
		for h in self.houses:
			x = h.power_interpolator(t) 
			self.pow_stats.append(x)
		for h in self.houses[:self.penetration]:
			try:
				x = h.generation_interpolator(t)
			except ValueError: 
				x = h.fill_missing_generation(t) 
			self.gen_stats.append(x) 

	def get_power_vec(self, t, **kwargs):
		self.__dict__.update(kwargs)
		self.Pvec[:] = 0.0 
		for i in range(self.n-1):
			self.Pvec[i] = self.houses[i].get_house_power(t) 
		surplus = sum(self.Pvec) 
		self.Pvec[-1] = -surplus
		return self.Pvec 

	def purge_batteries(self, **kwargs):
		self.__dict__.update(kwargs) 
		if self.batteries_on_grid:
			for h in self.houses:
				if h.has_battery:
					h.delete_battery() 
			self.batteries_on_grid = False 
		else:
			print("No batteries on grid.")



def bisect_into(v, x):
    lowerbound = 0
    upperbound = len(v) 
    while lowerbound < upperbound:
        midpoint = (lowerbound + upperbound)//2
        if x < v[midpoint]: 
        	upperbound = midpoint 
        else: 
        	lowerbound = midpoint + 1
    return lowerbound 



def select_random_element(v): 
	cumsum = np.cumsum(v) 
	total = cumsum[-1] 
	randomnum = rand()*total
	bisected_index = bisect_into(cumsum, randomnum) 
	return v[bisected_index] 



class Trajectory:    

	def __init__(self, **kwargs):
		self.__dict__.update(kwargs)
		self.n = 50 
		self.sigmas = []
		self.simplex_points = []
		self.maxpowers = []
		self.battery_level = []  
		self.ensemble_members = 0
		self.tot_pow = [] 
		self.tot_gen = [] 

	def add_tajectory_point(self, P, **kwargs):
		self.__dict__.update(kwargs) 
		source, sink, passive = continuoussourcesinkcounter(P)
		sigma = (source, sink, passive)
		self.sigmas.append(sigma) 
		self.maxpowers.append(np.max(np.abs(P))) 

	def add_battery_level(self, charge, **kwargs):
		self.__dict__.update(kwargs)
		self.battery_level.append(charge) 

	def convert_to_simplex_points(self, **kwargs):
		self.__dict__.update(kwargs)
		try:
			for sigma in self.sigmas: 
				s = sigma[0]*self.n
				d = sigma[1]*self.n 
				p = self.n - s - d 
				point = (s, d, p) 
				self.simplex_points.append(point) 
		except AttributeError:
			print("No data to convert")
			return  
		except IndexError: 
			print("No data to convert")
			return 

	def convert_a_point(self, point, **kwargs):
		self.__dict__.update(kwargs) 
		s = point[0]*self.n 
		d = point[1]*self.n 
		return (s, d, self.n-s-d)

	def update_mean(self, v, **kwargs):
		self.__dict__.update(kwargs)
		if self.ensemble_members == 0:
			self.mean_trajectory = [] 
			for triple in v: 
				sources, sinks, passives = triple[0], triple[1], triple[2]
				x = (sources, sinks, passives) 
				self.mean_trajectory.append(x) 
			self.ensemble_members += 1 
		else: 
			for i, triple in enumerate(v):
				z = self.mean_trajectory[i] 
				mean_sources, mean_sinks = z[0], z[1] 
				sources, sinks = triple[0], triple[1]
				mean_sources = (self.ensemble_members*mean_sources + sources)/(self.ensemble_members + 1.0) 
				mean_sinks = (self.ensemble_members*mean_sinks + sinks)/(self.ensemble_members + 1.0) 
				x = (mean_sources, mean_sinks, 1.0-mean_sources-mean_sinks) 
				self.mean_trajectory[i] = x 
			self.ensemble_members += 1 
		return self.mean_trajectory 

	def randomly_select_cascade_point(self, **kwargs):
		self.__dict__.update(kwargs)  
		return select_random_element(self.maxpowers)  





class TrajectoryPlotter: 

	n = 50

	def __init__(self, **kwargs):
		self.__dict__.update(kwargs)
		scale = self.n 
		self.figure, self.tax = ternary.figure(scale=scale-2)
		self.figure.set_size_inches(8, 8)
		self.figure.set_dpi(90)
		self.tax.boundary(linewidth=0.2) 

	def plot_trajectory(self, traj, alpha, lw, color, glinewidth, **kwargs):
		self.__dict__.update(kwargs)
		remapped_points = [] 
		for point in traj:
			ns = point[0]
			nd = point[1]
			ne = point[2] 
			x = (nd-1, ne, ns-1)
			remapped_points.append(x)
		self.tax.gridlines(linewidth=glinewidth, multiple=6, color="k")
		self.tax.plot(remapped_points, linewidth=lw, color=color, alpha=alpha)

	def plot_point(self, point, color="green", **kwargs): 
		self.__dict__.update(kwargs)
		ns = point[0]
		nd = point[1]
		ne = point[2] 
		x = [(nd-1, ne, ns-1)] 
		self.tax.scatter(x, color=color)


	def show_plot(self, title="simplex", **kwargs):
		self.__dict__.update() 
		self.tax.set_title(title)
		self.tax.get_axes().axis('off')
		self.tax.clear_matplotlib_ticks()
		self.tax.show()



























