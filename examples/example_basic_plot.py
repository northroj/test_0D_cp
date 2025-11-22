import h5py
import numpy as np
import matplotlib.pyplot as plt
with h5py.File("example_basic.h5","r") as f: 
    tally1 = f["tallies"]["user_tallies"]["test_tally_1"]
    energy_loss = np.array(tally1["counts"][:])
    energy_bins = np.array(tally1["energy_edges"][:])
    time_bins = np.array(tally1["time_edges"][:])
    x_bins = np.array(tally1["x_edges"][:])
    tally2 = f["tallies"]["user_tallies"]["test_tally_2"]
    average_particle_energy = np.array(tally2["counts"][:])
    average_p_e_time = np.array(tally2["time_edges"][:])
    """
    # standard tallies
    tallystandard = f["tallies"]["standard_tallies"]
    ion_temperature_time = np.array(tallystandard["ion_temperature_time"]["counts"][:])
    electron_temperature_time = np.array(tallystandard["electron_temperature_time"]["counts"][:])
    ion_density_time = np.array(tallystandard["ion_density_time"]["counts"][:])
    standard_time_bins = np.array(tallystandard["ion_density_time"]["time_edges"][:])
    """

"""
# standard tallies
plt.figure(1)
plt.plot(standard_time_bins[1:], ion_temperature_time[0,:,0], color="blue", label="ion")
plt.plot(standard_time_bins[1:], electron_temperature_time[0,:,0], color="red", label="electron")
plt.xlabel("Time [shk]")
plt.ylabel("Temperature [keV]")
plt.show()

plt.figure(3)
plt.plot(standard_time_bins[1:], ion_density_time[0,:,0], color="blue", label="d")
plt.plot(standard_time_bins[1:], ion_density_time[1,:,0], color="red", label="t")
plt.plot(standard_time_bins[1:], ion_density_time[2,:,0], color="green", label="a")
plt.xlabel("Time [shk]")
plt.ylabel("Density [g/cc]")
plt.yscale("log")
plt.show()
"""

# user tallies

plot_energy = np.sum(energy_loss[2,:,:], axis=0)
plot_time = np.sum(energy_loss[2,:,:], axis=1)
plt.figure(100)
plt.plot(energy_bins[1:], plot_energy)
plt.xlabel("Energy [keV]")
plt.ylabel("Energy loss [keV]")
plt.show()
plt.figure(101)
plt.plot(time_bins[1:], plot_time)
plt.xlabel("Time [shk]")
plt.ylabel("Energy loss [keV]")
plt.show()

kendra_data = np.loadtxt("avg_E_out")
kendra_data_2 = np.loadtxt("5gcc_avg_E_out")
plot_avg_particle_energy = average_particle_energy[0,:]
plt.figure(200)
plt.plot(average_p_e_time[1:], plot_avg_particle_energy, color="red", label="my code")
plt.plot(kendra_data[:,0], kendra_data[:,2], color="blue", label="kendra's data")
plt.plot(kendra_data_2[:,0], kendra_data_2[:,2], color="green", label="kendra's data updated")
plt.legend()
plt.xlabel("Time [shk]")
plt.ylabel("Average particle energy [keV]")
plt.xlim([0, 0.005])
#plt.ylim([])
plt.show()