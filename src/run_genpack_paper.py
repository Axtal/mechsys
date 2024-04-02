import subprocess

for sc in [0,1]: # for each sidecam
    subprocess.call(["./dem_genpack", "22", "1", "%d"%sc, "0"]) # num, aratio, sidecam, show
    subprocess.call(["./dem_genpack", "17", "2", "%d"%sc, "0"])
    subprocess.call(["./dem_genpack", "15", "3", "%d"%sc, "0"])
    subprocess.call(["./dem_genpack", "14", "4", "%d"%sc, "0"])
    subprocess.call(["./dem_genpack", "13", "5", "%d"%sc, "0"])
