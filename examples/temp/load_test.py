import flopy
pth = "/Users/langevin/langevin/dev/modflow6-testmodels.git/mf6/test014_NWTP3Low_dev"
sim = flopy.mf6.MFSimulation.load(sim_ws=pth)
new_pth = "./examples/temp/test014_NWTP3Low_dev"
sim.set_sim_path(new_pth)
sim.write_simulation()
print("done.")