import numpy as np 
import dill 
import matplotlib.pyplot as plt 
import glob, os
from powerclasses import * 
from matplotlib import colors 

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


# def compute_daily_profile_with_stats(n, month, penetration, selection1, selection2):
# 	G = MicroGrid() 
# 	G.n = n 
# 	G.penetration = penetration
# 	G.month = month
# 	G.make_houses() 
# 	G.assign_genprofiles()
# 	G.assign_powerprofiles() 

# 	#Move through the day 
# 	traj = Trajectory()
# 	traj.n = n  
# 	point1, point2 = 0, 0
# 	point1_powstats, point2_powstats, point1_genstats, point2_genstats = [], [], [], []  
# 	for t in tday_sample: 
# 		P = G.get_power_vec(t) 
# 		traj.add_tajectory_point(P)
# 		G.get_power_breakdown(t) 
# 		traj.tot_pow.append(G.pow_tot)
# 		traj.tot_gen.append(G.gen_tot)
# 		if t == selection1: 
# 			G.get_hists(t)
# 			point1 = traj.sigmas[-1]
# 			point1_genstats = np.array(G.gen_stats)  
# 			point1_powstats = np.array(G.pow_stats)  
# 		if t == selection2:
# 			G.get_hists(t)
# 			point2 = traj.sigmas[-1] 
# 			point2_genstats = np.array(G.gen_stats)  
# 			point2_powstats = np.array(G.pow_stats)  
# 	traj.convert_to_simplex_points() 
# 	point1 = traj.convert_a_point(point1) 
# 	point2 = traj.convert_a_point(point2) 
# 	return traj, point1_genstats, point1_powstats, point2_genstats, point2_powstats, point1, point2 

# tr, p1g, p1p, p2g, p2p, p1, p2 = compute_daily_profile_with_stats(20, 1, 19, tday_sample[2], tday_sample[24])
# plot = TrajectoryPlotter(n=20)
# plot.plot_trajectory(tr.simplex_points, 0.7, 1.5, "k", 0.05)
# plot.plot_point(p1, color="green")
# plot.plot_point(p2, color="red") 
# plot.show_plot() 
# plt.plot(tday_sample, tr.tot_pow)
# plt.plot(tday_sample, tr.tot_gen)
# plt.axvline(tday_sample[2])
# plt.axvline(tday_sample[24]) 
# plt.show() 
# plt.clf() 
# xinds = np.arange(0,19,1)
# print(p1p)
# plt.bar(xinds, p1p, color=martaRed)
# plt.bar(xinds, p1g, bottom=p1p, color=martaBlue)
# plt.show() 
# plt.clf()
# plt.bar(xinds, p2p, color=martaRed)
# plt.bar(xinds, p2g, bottom=p2p, color=martaBlue)
# plt.show() 
# np.savetxt("data/p1p.txt", p1p)
# np.savetxt("data/p1g.txt", p1g)
# np.savetxt("data/p2p.txt", p2p)
# np.savetxt("data/p2g.txt", p2g) 


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


def read_in_ensemble_and_make_means(directory, n): 

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

	plot.show_plot() 


# do_ensemble_trajectories(50, 50, 1, 24, 49)

# tr1, tr2, tr3, tr4 = compute_all_weekly_profiles(50, 1, 24, 49) 
# plot = TrajectoryPlotter(n=50)
# plot.plot_trajectory(tr1.simplex_points, 0.2, 3.0, "k", 0.05)
# plot.plot_trajectory(tr2.simplex_points, 0.2, 3.0, "k", 0.05)
# plot.plot_trajectory(tr3.simplex_points, 0.2, 3.0, "k", 0.05)
# plot.plot_trajectory(tr4.simplex_points, 0.2, 3.0, "k", 0.05)
# plot.show_plot() 
# # plt.plot(traj.maxpowers)
# plt.show()


read_in_ensemble_and_make_means("trajdata2/7/fullpen_withbat/", 50)



# with open("data/trajdata/10/fullpen_nobat/8.pkl", "rb") as f:
# 	traj = dill.load(f, encoding='latin1') 
# print(type(traj))
# print(traj.n)
# for x in traj.maxpowers:
# 	print(x)
# plot = TrajectoryPlotter(n=traj.n)
# plot.plot_trajectory(traj.simplex_points, 0.5, 1.5, "k", 0.05)
# plot.show_plot() 




# with open("data/trajdata/1/fullpen_nobat/0.pkl", "rb") as f:
# 	traj = dill.load(f) 
# print(type(traj))
# print(traj.n)
# for x in traj.maxpowers:
# 	print(x)
# plot = TrajectoryPlotter(n=100)
# plot.plot_trajectory(traj.simplex_points, 0.6, 3.0, "k")
# plot.show_plot() 
# plt.plot(traj.maxpowers)
# plt.show()





# G = MicroGrid(n=2, month=1, penetration=1) 
# G.make_houses() 
# G.assign_genprofiles()
# G.assign_powerprofiles() 
# G.assign_batteries() 
# d, p, s, ptot = [], [], [], []
# for t in tweek_sample:
# 	selected_house = G.houses[0]
# 	P = selected_house.get_house_power(t)
# 	ptot.append(P)
# 	demand, production, storage = selected_house.get_states(t)
# 	d.append(-demand)
# 	p.append(production)
# 	s.append(storage)

# plt.plot(d, color="k")
# plt.plot(p, color="r")
# plt.plot(s, "o-", ms=3, color="b")
# plt.plot(ptot, "k--", lw=3.0, alpha=0.5)
# plt.show()



# tr1, tr2 = compute_both_weekly_profiles(10, 1, 4)
# # plot = TrajectoryPlotter(n=10)
# # plot.plot_trajectory(tr1.simplex_points, 0.6, 3.0, "k")
# # # plot.plot_trajectory(tr2.simplex_points, 0.6, 3.0, "b")
# # plot.show_plot() 

# print("none batteries:")
# for point in tr1.battery_level:
# 	print(point)
# print("batteries:")
# for point in tr2.battery_level:
# 	print(point)

# plt.plot(tr1.battery_level)
# plt.plot(tr2.battery_level, color="r")
# plt.show() 




# traj = compute_weekly_profile(6, 1, 6)
# for x in traj.simplex_points:
# 	print(x)
# print("Powers: ")
# for x in traj.maxpowers:
# 	print(x)  

# plot = TrajectoryPlotter(n=6)
# plot.plot_trajectory(traj.simplex_points, 0.6, 3.0)
# plot.show_plot() 

# with open("data/trajectories/test.pkl", "wb") as f:
# 	dill.dump(traj, f) 
# with open("data/trajectories/test.pkl", "rb") as f:
# 	traj = dill.load(f) 
# print(type(traj))
# for x in traj.maxpowers:
# 	print(x)
# plot = TrajectoryPlotter(n=6)
# plot.plot_trajectory(traj.simplex_points, 0.6, 3.0)
# plot.show_plot() 






















