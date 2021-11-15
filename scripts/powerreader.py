import dask.dataframe as dd
import glob, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates 
import matplotlib
from matplotlib import colors
import datetime as dt
from datetime import timedelta 
import glob, os
from numpy.random import randint 

#Set up plot style
font = {'size'   : 14}
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


def gettimesandpower(date, houseid, fname):
	"""
	For a given houseid and time (or time window), the function
	returns two vectors: the times and power levels, respectively.
	"""
	df = dd.read_csv(fname) 

	#Filter for id
	x = df.loc[df["LCLid"] == houseid]
	x = x.loc[x["DateTime"].str.contains(date), ["DateTime", "KWH/hh (per half hour) "]]
	y = x.compute() 

	#Format
	y["DateTime"] = y["DateTime"].map(lambda d: dt.datetime.strptime(d[:-8],'%Y-%m-%d %H:%M:%S'))
	y["KWH/hh (per half hour) "] = y["KWH/hh (per half hour) "].map(lambda d: float(d))
	y = np.array(y) 
	return y[:,0], y[:,1] 


def getpowerattimes(date, time, houseid, fname):
	"""
	For a given houseid and set of dates, the function returns
	the power levels at a given time on those dates.
	"""
	df = dd.read_csv(fname) 

	#Filter for id
	x = df.loc[df["LCLid"] == houseid]
	x = df.loc[df["DateTime"].str.contains(date)]
	x = x.loc[x["DateTime"].str.contains(time), ["KWH/hh (per half hour) "]]
	y = x.compute() 
	y["KWH/hh (per half hour) "] = y["KWH/hh (per half hour) "].map(lambda d: float(d))
	y = np.array(y)
	return y


def PVgettimesandpower(date, houseid, fname):
	"""
	For a given houseid and time (or time window), the function
	returns two vectors: the times and power levels, respectively.
	"""
	df = dd.read_csv(fname)
	x = df.loc[df["Substation"] ==  houseid]
	x = x.loc[x["datetime"].str.contains(date), ["datetime", "P_GEN"]]
	y = x.compute() 
	y["datetime"] = y["datetime"].map(lambda d: dt.datetime.strptime(d,'%Y-%m-%d %H:%M:%S'))
	y["P_GEN"] = y["P_GEN"].map(lambda d: float(d))
	y = np.array(y) 
	return y[:,0], y[:,1] 


def trimandshift(t, p):
	"""
	Takes the vectors given by the above function, 
	and chops out a day and samples at half hourly intervals.
	"""
	t3, p3 = [], [] 
	t2 = t[132+4*144:276+4*144]
	p2 = p[132+4*144:276+4*144]
	for i in range(len(t2)): 
		if (i%3) == 0: 
			t3.append(t2[i])
			p3.append(p2[i])
	return t3, p3 


def getPVvecs(fname):
	"""
	Generates an ensemble of day long PV activities, sampled  3 different 
	days for each complete pv data set 
	"""
	datmat = np.zeros((18,48))
	df = dd.read_csv(fname)
	i = 0 
	for unique_value in df.Substation.unique():
		ttemp, ptemp = PVgettimesandpower("2014-06", unique_value, fname)
		t, p = trimandshift(ttemp, ptemp) 
		datmat[i,:] = np.array(p) 
		i += 1
		ttemp, ptemp = PVgettimesandpower("2014-07", unique_value, fname)
		t, p = trimandshift(ttemp, ptemp) 
		datmat[i,:] = np.array(p) 
		i += 1
		ttemp, ptemp = PVgettimesandpower("2014-08", unique_value, fname)
		t, p = trimandshift(ttemp, ptemp) 
		datmat[i,:] = np.array(p) 
		i += 1
	return datmat 


def sourcesinkcounter(Pvec, tol=0.1):
	ns = 0 
	nd = 0
	for i in range(len(Pvec)):
		if Pvec[i] > tol: 
			ns += 1 
		if Pvec[i] < -tol:
			nd += 1 
	ne = len(Pvec) - ns - nd  
	return ns, nd, ne 


def continuoussourcesinkcounter(Pvec):
	largestsource = np.max(Pvec)
	n = len(Pvec)
	largestsink = np.abs(np.min(Pvec)) 
	sourceterms = list(filter(lambda x: x>0.0, Pvec))
	sinkterms = list(filter(lambda x: x<0.0, Pvec))
	sigma_s = np.sum(sourceterms)/(n*largestsource)
	sigma_d = np.sum(np.abs(sinkterms))/(n*largestsink) 
	sigma_p = 1.0 - sigma_s - sigma_d 
	return sigma_s, sigma_d, sigma_p 
	

def import_data():
	raw_data_df = pd.read_csv("powerdata/data/Power-Networks-LCL-June2015(withAcornGps)v2_1.csv", header=0)
	return raw_data_df


def strip_to_datevpower(res):
	res["date"] = pd.to_datetime(res["DateTime"])
	data = res.loc[:, ["KWH/hh (per half hour) "]] 
	data = data.set_index(res.date)
	data["KWH/hh (per half hour) "] = pd.to_numeric(data["KWH/hh (per half hour) "], downcast="float", errors="coerce")
	return data 

def strip_to_datevidandpower(res):
	res["date"] = pd.to_datetime(res["DateTime"])
	data = res.loc[:, ["LCLid", "KWH/hh (per half hour) "]]
	data = data.set_index(res.date)
	data["KWH/hh (per half hour) "] = pd.to_numeric(data["KWH/hh (per half hour) "], downcast="float", errors="coerce")
	return data 


def resampler(dat, period):
	if period == "weekly":
		newdat = dat.resample("W").sum()
	else:
		newdat = dat.resample("H").sum() 
	return newdat 


def get_unique_house_power(raw_df, identifier):
		results = raw_df.loc[raw_df["LCLid"] ==  identifier] 
		data = strip_to_datevpower(results) 


def get_house_power_sum(raw_df, ax):
	i = 0 
	for val in raw_df.LCLid.unique(): 
		if i == 3:
			print(val)
			results = raw_df.loc[raw_df["LCLid"] ==  val] 
			data = strip_to_datevpower(results) 
			thinned_data = resampler(data, "weekly")
			# print("Stuff: ", thinned_data.loc["2011-12-11", "KWH/hh (per half hour) "])
			ax.plot(thinned_data, color=martaRed, lw=2.0, alpha=0.9)
			return
		else: 
			i += 1


def f2(raw_df):
	"""
	Returns dataframe of (date, power, ID) but thinned to weekly. 
	Different IDs appended head to tail.
	"""
	df = strip_to_datevidandpower(raw_df) 
	for (i, val) in enumerate(df.LCLid.unique()):
		if i == 0: 
			selection = df.loc[df["LCLid"] == val] 
			thinned  =  resampler(selection, "d")
			thinned["LCLid"] = val 
		else:
			selection = df.loc[df["LCLid"] == val]
			temp = resampler(selection, "d")
			temp["LCLid"] = val
			thinned = pd.concat([thinned,temp], ignore_index=False)
	return thinned


def meanyears(thinned_df): 
	dates, stds, means, ns = [], [], [], []
	for datevals in thinned_df.index.unique():
		pows = thinned_df.loc[datevals, "KWH/hh (per half hour) "].to_numpy()
		means.append(np.mean(pows))
		stds.append(np.std(pows))
		dates.append(datevals)
		ns.append(len(pows))
	return dates, stds, means, ns



def expand_mean_and_std(n_prev, mean_prev, std_prev, newdat):
	n_new = len(newdat) 
	mean_new = np.mean(newdat) 
	std_new = np.std(newdat) 
	mean_combined = (n_prev*mean_prev + n_new*mean_new)/(n_prev+n_new)
	d1 = mean_prev - mean_combined 
	d2 = mean_new - mean_combined
	v1 = std_prev**2 
	v2 = std_new**2  
	v_combined = (v1/n_prev) + (v2/n_new)
	return mean_combined, np.sqrt(v_combined)



def wholeset_meanyears():

	#Crawl through the files
	search_dir = "powerdata/data/"
	files = filter(os.path.isfile, glob.glob(search_dir + "*.csv"))
	new_df = [] 
	for i, f in enumerate(files):
		print(i, f)
		if i==0:
			raw_data_df = pd.read_csv(f, header=0)
			thinned_df = f2(raw_data_df) 
			d, s, m, n = meanyears(thinned_df)
			dat = np.vstack((m,s,n)).transpose()
			new_df = pd.DataFrame(dat, index=d, columns=["means","stds","n"]) 
			new_df.index= pd.to_datetime(new_df.index) 
		else:
			raw_data_df = pd.read_csv(f, header=0)
			thinned_df = f2(raw_data_df) 
			for dateval in thinned_df.index.unique():
				if dateval in new_df.index:

					#Get the values of mean, std and n from df 
					mean_prev, std_prev, n_prev = new_df.loc[dateval].to_numpy()

					try: 
						#Put the new data into the expanded mean function
						pows = thinned_df.loc[dateval, "KWH/hh (per half hour) "]
						pows = pows.to_numpy()
						mean_comb, std_comb = expand_mean_and_std(n_prev, mean_prev, std_prev, pows)
						n_comb = n_prev + len(pows)

						#Update the df values of mean, std and n 
						new_df.loc[dateval] = [mean_comb, std_comb, n_comb]
					except AttributeError:
						pass 

				else:
					try: 
						pows = thinned_df.loc[dateval, "KWH/hh (per half hour) "]
						pows = pows.to_numpy() 
						mean = np.mean(pows) 
						std = np.std(pows) 
						print("Inserting...")
						new_df.loc[dateval] = [mean, std, len(pows)] 
					except AttributeError:
						pass 	
	return new_df


def wholeset_meanyears_PV():

	raw_data_df = pd.read_csv("powerdata/PV/HourlyData/CustomerEndpoints.csv")

	#Strip to date vs id and power
	raw_data_df.datetime = pd.to_datetime(raw_data_df.datetime) 
	data = raw_data_df.loc[:, ["Substation", "P_GEN_MAX", "P_GEN_MIN"]]
	data = data.set_index(raw_data_df.datetime)
	data.P_GEN_MAX = pd.to_numeric(data.P_GEN_MAX, downcast="float", errors="coerce")
	data.P_GEN_MIN = pd.to_numeric(data.P_GEN_MIN, downcast="float", errors="coerce")

	#Re-sample to weekly
	for (i, val) in enumerate(data.Substation.unique()): 
		if i == 0: 
			selection = data.loc[data.Substation == val] 
			thinned_df =  resampler(selection, "h")
			thinned_df.Substation= val 
		else:
			selection = data.loc[data.Substation == val]
			temp = resampler(selection, "h")
			temp.Substation = val
			thinned_df = pd.concat([thinned_df,temp], ignore_index=False) 

	dates, stds, means, ns = [], [], [], []
	for datevals in thinned_df.index.unique():
		try:
			pows1 = thinned_df.loc[datevals, "P_GEN_MAX"].to_numpy()
			pows2 = thinned_df.loc[datevals, "P_GEN_MIN"].to_numpy()
			pows = (pows1 + pows2)/2.0 
			means.append(np.mean(pows))
			stds.append(np.std(pows))
			dates.append(datevals)
			ns.append(len(pows))
		except AttributeError:
			pows1 = thinned_df.loc[datevals, "P_GEN_MAX"]
			pows2 = thinned_df.loc[datevals, "P_GEN_MIN"]
			print(pows1, pows2)
			pows = np.mean((pows1,pows2))
			print(pows)
			std01 = np.std((pows1,pows2))
			print(std01)
			print(datevals)
			means.append(pows)
			stds.append(std01)
			dates.append(datevals)
			ns.append(1)
	#Sort into time order, should be mostly there already
	z = np.array(sorted(zip(dates,means,stds,ns), key=lambda pair: pair[0])) 
	dates, means, stds, ns = z[:,0], z[:,1], z[:,2], z[:,3] 

	return dates, stds, means, ns




monthlist = ["January","February","March","April","May","June","July",
	"August","September","October","November","December"]
	

def make_random_week_profiles(monthchoice):

	#Choose random file
	search_dir = "powerdata/data/"
	files = filter(os.path.isfile, glob.glob(search_dir + "*.csv"))
	filelist = list(files)
	numfiles = len(filelist) 
	filechoice = randint(0,numfiles)

	#Open file
	f = filelist[filechoice]
	raw_df = pd.read_csv(f, header=0)
	df = strip_to_datevidandpower(raw_df)

	#Choose a random house 
	houeselist = list(df.LCLid.unique())
	numhouses = len(houeselist)
	housechoice = randint(0,numhouses) 
	house = houeselist[housechoice] 

	#Open data for that house
	selection = df.loc[df["LCLid"] == house] 
	selection.index = pd.to_datetime(selection.index) 

	#Choose the desired month
	month = selection.loc[selection.index.to_series().dt.month==monthchoice] 

	#Find days of week for each day in month 
	month["day_of_week"] = month.index.to_series().dt.dayofweek
	month = month.sort_index() 

	#Check a full week exists 
	if len(month.day_of_week.unique()) == 7:

		#Set artificial dates to each day 
		#so that a mean time series can be built
		indexaslist = month.index.tolist() 
		for (p, val) in enumerate(month.index):
			for l in range(0,7):
				if month.iloc[p].day_of_week == l:
					indexaslist[p] = indexaslist[p].replace(year=1970,month=1,day=l+1)
		month.index = indexaslist

		#Now construct the means for each time
		dates, stds, means = [], [], []
		for l in range(0,7):
			pows = month.loc[month["day_of_week"]==l]
			for t in pows.index.unique():
				pows2 = pows.loc[pows.index==t, "KWH/hh (per half hour) "]
				try:
					pows2 = pows2.to_numpy()
					means.append(np.mean(pows2))
					stds.append(np.std(pows2))
				except AttributeError:
					pows2 = pows.loc[datevals, "KWH/hh (per half hour) "]
					means.append(pows2)
					stds.append(0.0)
				dates.append(t)

		#Sort into time order, should be mostly there already
		z = np.array(sorted(zip(dates,means,stds), key=lambda pair: pair[0])) 
		dates, means, stds = z[:,0], z[:,1], z[:,2] 

		#Make a seconds axis
		dates2 = pd.Series(dates)
		idx = pd.to_timedelta(dates2)
		secs = np.array(idx.dt.total_seconds())

		#Output weekly time series 
		outdf = pd.DataFrame({"secs":secs,"means":means,"stds":stds}, index=dates) 
		return outdf, True
	else:
		return 0.0, False 


def make_random_week_profiles_PV(monthchoice):

	#Open file
	f = "powerdata/PV/HourlyData/CustomerEndpoints.csv"
	raw_df = pd.read_csv(f, header=0)

	#Strip data 
	raw_df.datetime = pd.to_datetime(raw_df.datetime) 
	df = raw_df.loc[:, ["Substation", "P_GEN_MAX", "P_GEN_MIN"]]
	df = raw_df.set_index(raw_df.datetime)
	df.P_GEN_MAX = pd.to_numeric(df.P_GEN_MAX, downcast="float", errors="coerce")
	df.P_GEN_MIN = pd.to_numeric(df.P_GEN_MIN, downcast="float", errors="coerce")

	#Choose a random PV panel 
	pvlist = list(df.Substation.unique())
	numpv = len(pvlist)
	pvchoice = randint(0,numpv) 
	pv = pvlist[pvchoice] 

	#Open data for that panel
	selection = df.loc[df["Substation"] == pv] 
	selection.index = pd.to_datetime(selection.index) 

	#Choose the desired month
	month = selection.loc[selection.index.to_series().dt.month==monthchoice] 

	#Find days of week for each day in month 
	month["day_of_week"] = month.index.to_series().dt.dayofweek
	month = month.sort_index() 

	#Check a full week exists 
	if len(month.day_of_week.unique()) == 7:

		#Set artificial dates to each day 
		#so that a mean time series can be built
		indexaslist = month.index.tolist() 
		for (p, val) in enumerate(month.index):
			for l in range(0,7):
				if month.iloc[p].day_of_week == l:
					indexaslist[p] = indexaslist[p].replace(year=1970,month=1,day=l+1)
		month.index = indexaslist

		#Now construct the means for each time
		dates, stds, means = [], [], []
		for l in range(0,7):
			pows = month.loc[month["day_of_week"]==l]
			for t in pows.index.unique():
				pows2 = pows.loc[pows.index==t, "P_GEN_MAX"].to_numpy()
				pows3 = pows.loc[pows.index==t, "P_GEN_MIN"].to_numpy()
				pows2 = (pows2 + pows3)/2.0 
				try:
					means.append(np.mean(pows2))
					stds.append(np.std(pows2))
				except AttributeError:
					means.append(pows2)
					stds.append(0.0)
				dates.append(t) 

		#Sort into time order, should be mostly there already
		z = np.array(sorted(zip(dates,means,stds), key=lambda pair: pair[0])) 
		dates, means, stds = z[:,0], z[:,1], z[:,2] 

		#Make a seconds axis
		dates2 = pd.Series(dates)
		idx = pd.to_timedelta(dates2)
		secs = np.array(idx.dt.total_seconds())

		#Output weekly time series 
		outdf = pd.DataFrame({"secs":secs,"means":means,"stds":stds}, index=dates) 
		return outdf, True
	else:
		return 0.0, False 
			













