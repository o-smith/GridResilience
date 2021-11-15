import numpy as np 
import dill 
import glob, os
from powerclasses import * 


t = np.linspace(0,604800-1800,336)[:-24]
tday_sample = t[48:96] 
tweek_sample = t[48:] 

#Function to compute a daily profile
def compute_daily_profile(n, month, penetration): 
	G = MicroGrid(n=n, month=month, penetration=penetration) 
	G.make_houses() 
	G.assign_genprofiles()
	G.assign_powerprofiles() 

	#Move through the day 
	traj = Trajectory(n=n) 
	for t in tday_sample: 
		P = G.get_power_vec(t) 
		traj.add_tajectory_point(P)
	traj.convert_to_simplex_points() 
	return traj 

#Function to compute a weekly profile 
def compute_weekly_profile(n, month, penetration, batteries=False): 
	G = MicroGrid(n=n, month=month, penetration=penetration) 
	G.make_houses() 
	G.assign_genprofiles()
	G.assign_powerprofiles() 

	#Move through the day 
	traj = Trajectory(n=n) 
	for t in tweek_sample: 
		P = G.get_power_vec(t) 
		traj.add_tajectory_point(P)
	traj.convert_to_simplex_points() 
	return traj 


def compute_all_weekly_profiles(n, month, pen1, pen2): 
	G = MicroGrid()
	G.n = n 
	G.month = month 
	G.penetration = pen1 
	G.make_houses() 
	G.assign_genprofiles()
	G.assign_powerprofiles() 
	print("Time-series assigned")

	#Move through the day 
	print("Standard case ...")
	traj1 = Trajectory()
	traj1.n = n 
	for t in tweek_sample: 
		P = G.get_power_vec(t) 
		mean_charge = G.get_av_battery_charge() 
		traj1.add_tajectory_point(P)
		traj1.add_battery_level(mean_charge)
	traj1.convert_to_simplex_points()

	#Now do the same again but with batteries on 
	print("Battery case ...")
	G.assign_batteries() 
	traj2 = Trajectory()
	traj2.n = n  
	for t in tweek_sample: 
		P = G.get_power_vec(t) 
		mean_charge = G.get_av_battery_charge() 
		traj2.add_tajectory_point(P)
		traj2.add_battery_level(mean_charge)
	traj2.convert_to_simplex_points()

	#Now repeat the above for the second level of penetration 
	print("Starting second penetration level ...") 
	G.purge_batteries() 
	G.penetration = pen2 
	G.assign_genprofiles()
	print("Standard case ...")
	traj3 = Trajectory()
	traj3.n = n  
	for t in tweek_sample: 
		P = G.get_power_vec(t) 
		mean_charge = G.get_av_battery_charge() 
		traj3.add_tajectory_point(P)
		traj3.add_battery_level(mean_charge)
	traj3.convert_to_simplex_points() 

	G.assign_batteries()
	traj4 = Trajectory()
	traj4.n = n 
	for t in tweek_sample: 
		P = G.get_power_vec(t) 
		mean_charge = G.get_av_battery_charge() 
		traj4.add_tajectory_point(P)
		traj4.add_battery_level(mean_charge)
	traj4.convert_to_simplex_points()
	del G

	#Outputs:
	#    traj1: pen level 1, no batteries 
	#    traj2: pen level 1, with batteries
	#    traj3: pen level 2, no batteries 
	#    traj4: pen level 2, with batteries 
	return traj1, traj2, traj3, traj4  


def do_ensemble_trajectories(ensemble_size, n, month, pen1, pen2):
	for z in range(ensemble_size):
		print("Ensemble number: ", z)
		tr1, tr2, tr3, tr4 = compute_all_weekly_profiles(n, month, pen1, pen2) 
		root_dir = "data/trajdata2/"
		month_dir = str(month) + "/" 
		fname_tr1 = root_dir + month_dir + "halfpen_nobat/" + str(z) + ".pkl"
		fname_tr2 = root_dir + month_dir + "halfpen_withbat/" + str(z) + ".pkl"
		fname_tr3 = root_dir + month_dir + "fullpen_nobat/" + str(z) + ".pkl"
		fname_tr4 = root_dir + month_dir + "fullpen_withbat/" + str(z) + ".pkl"
		with open(fname_tr1, "wb") as f:
			dill.dump(tr1, f)
		with open(fname_tr2, "wb") as f:
			dill.dump(tr2, f)
		with open(fname_tr3, "wb") as f:
			dill.dump(tr3, f)
		with open(fname_tr4, "wb") as f:
			dill.dump(tr4, f)
		del tr1, tr2, tr3, tr4 
	return


def read_in_ensemble_and_make_means(directory, title, n): 

	plot = TrajectoryPlotter(n=n)

	#Crawl through the files
	files = filter(os.path.isfile, glob.glob(directory + "*.pkl")) 
	mean_traj = Trajectory() 
	for i, f in enumerate(files): 
		print(i) 

		#Un-pickle the file 
		with open(f, "rb") as opened_file:
			traj = dill.load(opened_file, encoding="latin1") 
			traj.convert_to_simplex_points()
		if i==0:
			plot.plot_trajectory(traj.simplex_points, 0.02, 0.9, "#404040", 0.03) 
		else: 
			plot.plot_trajectory(traj.simplex_points, 0.02, 0.9, martaRed, 0.0) 
		mean_traj.update_mean(traj.sigmas) 
		del(traj)
	mean_traj.sigmas = mean_traj.mean_trajectory  
	mean_traj.convert_to_simplex_points() 
	plot.plot_trajectory(mean_traj.simplex_points, 0.8, 2.8, martaRed, 0.0) 

	plot.show_plot(title) 

try:
	#Computing ensembles of random micro-grids from 
	#the power usage data for half and full uptake, with and
	#without batteries in January, April, July and October
	do_ensemble_trajectories(50, 50, 1, 24, 49) 
	do_ensemble_trajectories(50, 50, 4, 24, 49) 
	do_ensemble_trajectories(50, 50, 7, 24, 49) 
	do_ensemble_trajectories(50, 50, 10, 24, 49) 
except FileNotFoundError: 
	#PV and/or generation data has not been downloaded, 
	#use pre-computed trajectories 
	print("Reverting to precomputed trajectories...")

print("Plotting trajectories...")
read_in_ensemble_and_make_means("trajdata/1/halfpen_nobat/", "Winter: 50% penetration, no batteries", 50)
read_in_ensemble_and_make_means("trajdata/1/fullpen_nobat/", "Winter: 100% penetration, no batteries", 50)
read_in_ensemble_and_make_means("trajdata/1/halfpen_withbat/", "Winter: 100% penetration, with batteries", 50)
read_in_ensemble_and_make_means("trajdata/7/halfpen_nobat/", "Summer: 50% penetration, no batteries", 50)
read_in_ensemble_and_make_means("trajdata/7/fullpen_nobat/", "Summer: 100% penetration, no batteries", 50)
read_in_ensemble_and_make_means("trajdata/7/halfpen_withbat/", "Summer: 100% penetration, with batteries", 50)
























