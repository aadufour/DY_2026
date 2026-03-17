import pylhe
import numpy as np
import matplotlib.pyplot as plt

lhe_file_1 = "/Users/albertodufour/MG5/mg5amcnlo/DY_2/myLHE/unweighted_events.lhe"
lhe_file_2 = "/Users/albertodufour/MG5/mg5amcnlo/all_contributing/myLHE/unweighted_events.lhe"


def mll(p1, p2):
    p = np.array(p1) + np.array(p2)
    return np.sqrt(p[3]**2 - sum(p[l]**2 for l in range(3)))

def rap(p1, p2):
    p = np.array(p1) + np.array(p2)
    E  = p[3]
    pz = p[2]
    return 0.5 * np.log((E + pz) / (E - pz))

def cstar(p1, p2):
    # p1 = e-, p2 = e+
    p1 = np.array(p1)
    p2 = np.array(p2)
    
    # dilepton system
    p = p1 + p2
    E  = p[3]
    pz = p[2]
    mass = mll(p1, p2)  # or compute inline
    
    # boost the e- to the dilepton rest frame
    # boost vector along z
    beta  = pz / E
    gamma = E / mass
    
    # boosted e- momentum components
    pz1_boosted = gamma * (p1[2] - beta * p1[3])
    E1_boosted  = gamma * (p1[3] - beta * p1[2])
    
    # magnitude of 3-momentum of e- in rest frame
    p1_mag = np.sqrt(p1[0]**2 + p1[1]**2 + pz1_boosted**2)
    
    # direction of dilepton motion in lab = z axis
    # cos theta* is angle between e- and z in rest frame
    return pz1_boosted / p1_mag
    



if __name__ == "__main__":
    inv_mass_1 = []
    rapidity_1 = []
    cosstar_1 = []

    inv_mass_2 = []
    rapidity_2 = []
    cosstar_2 = []

    i=0
    for event_1, event_2 in zip(pylhe.read_lhe_with_attributes(lhe_file_1), pylhe.read_lhe_with_attributes(lhe_file_2)):
        i += 1
        if not (i)%1000:
            print(f"\n--- Event {i} ---")
        final_state_1 = [particle for particle in event_1.particles if particle.status == 1.0]
        pmu = [[final_state_1[k].px, final_state_1[k].py, final_state_1[k].pz, final_state_1[k].e] for k in range(2)]
        inv_mass_1.append(mll(pmu[0], pmu[1]))
        rapidity_1.append(rap(pmu[0], pmu[1]))
        cosstar_1.append(cstar(pmu[0], pmu[1]))

        final_state_2 = [particle for particle in event_2.particles if particle.status == 1.0]
        pmu = [[final_state_2[k].px, final_state_2[k].py, final_state_2[k].pz, final_state_2[k].e] for k in range(2)]
        inv_mass_2.append(mll(pmu[0], pmu[1]))
        rapidity_2.append(rap(pmu[0], pmu[1]))
        cosstar_2.append(cstar(pmu[0], pmu[1]))
    










    plt.xlim(50, 130)
    plt.hist(inv_mass_1, bins=2000, alpha=0.3, color="red", label="sm")
    plt.hist(inv_mass_2, bins=2000, alpha=0.3, color="blue", label="eft")

    
    plt.xlabel("m_ll")
    plt.legend()
    plt.savefig('/Users/albertodufour/code/DY2026/test/mll_dist_smeftsim.png')
   
    plt.clf()

    plt.hist(rapidity_1, bins=100, alpha=0.3, color="red", label="sm")
    plt.hist(rapidity_2, bins=100, alpha=0.3, color="blue", label="eft")
    plt.xlabel("y")
    plt.legend()
    plt.savefig('/Users/albertodufour/code/DY2026/test/rap_dist_smeftsim.png')
    plt.clf()

    plt.hist(cosstar_1, bins=100, alpha=0.3, color="red", label="sm")
    plt.hist(cosstar_2, bins=100, alpha=0.3,  color="blue", label="eft")
    plt.xlabel("c*")
    plt.legend()
    plt.savefig('/Users/albertodufour/code/DY2026/test/cstar_dist_smeftsim.png')

        