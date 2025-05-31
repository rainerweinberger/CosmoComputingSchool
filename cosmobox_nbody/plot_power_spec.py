import numpy as np
import matplotlib.pyplot as plt
import os

fout = "powerspec_z0.png"
RunDir = "./output/"
Tag = "DM-L50-N32"
inpsec_filename = "./output/inputspec_snapshot.txt"

BoxSize = 50.0
Num = 2

exts = f"{Num:04d}" if Num >= 1000 else f"{Num:03d}"

# Prepare plot
plt.figure(figsize=(10, 7.5))
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$k\ [h\ \mathrm{Mpc}^{-1}]$")
plt.ylabel(r"$\Delta^2(k)$")
plt.xlim(2 * np.pi / BoxSize * 0.7, 2 * np.pi / BoxSize * 2000)
plt.ylim(0.1, 8.0e3)

Kall = []
Dall = []

FoldFac = 16
MinModeCount = 8
TargetBinNummer = 50

fname = os.path.join(RunDir, f"powerspecs/powerspec_{exts}.txt")

with open(fname, 'r') as file:
    for piece in range(3):
        Time = float(file.readline())
        Bins = int(file.readline())
        BoxSize = float(file.readline())
        Ngrid = float(file.readline())
        Dgrowth = float(file.readline())
        
        da = np.loadtxt(file, max_rows=Bins).T
        K, Delta2, ModePow, ModeCount, Shot = da

        if piece == 0:
            with open(inpsec_filename, 'r') as f2:
                z, dplus = f2.readline().strip().split()
                dat = np.loadtxt(f2, max_rows=514).T
            k_lin, D2_lin = dat[0], dat[1] / Dgrowth**2
            plt.plot(k_lin, D2_lin)

        plt.plot(K, Shot, linestyle='--', linewidth=3.0)
        Delta2 -= Shot

        kmin = 2 * np.pi / BoxSize * FoldFac**piece
        kmax = 2 * np.pi / BoxSize * Ngrid / 2.0 / 4 * FoldFac**piece
        if piece > 0:
            kmin = 2 * np.pi / BoxSize * Ngrid / 2.0 / 4 * FoldFac**(piece - 1)

        MinDlogK = (np.log10(K.max()) - np.log10(K.min())) / TargetBinNummer

        istart = 0
        k_list, Delta2_list = [], []

        while istart < Bins:
            ind = slice(istart, istart + 1)
            count = ModeCount[ind].sum()
            deltak = np.log10(K[ind].max()) - np.log10(K[ind].min())

            if deltak >= MinDlogK:
                weights = ModeCount[ind]
                kk = np.exp(np.sum(np.log(K[ind]) * weights) / weights.sum())
                d2 = np.sum(Delta2[ind] * weights) / weights.sum()
                k_list.append(kk)
                Delta2_list.append(d2)
                istart += 1
            else:
                istart += 1

        mask = (np.array(k_list) >= kmin) & (np.array(k_list) <= kmax)
        Kall.extend(np.array(k_list)[mask])
        Dall.extend(np.array(Delta2_list)[mask])

Kall = np.array(Kall)
Dall = np.array(Dall)

Sall = np.interp(np.log(Kall), np.log(K), np.log(Shot))
Sall = np.exp(Sall)

mask = Dall > Sall / 2
plt.plot(Kall[mask], Dall[mask], color='black', linewidth=4)

plt.text(0.22, 0.76, Tag, transform=plt.gca().transAxes, color='green')
plt.text(0.22, 0.85, "power spectrum, z = 0", transform=plt.gca().transAxes)

plt.savefig(fout, format='png', dpi=200)
plt.close()
