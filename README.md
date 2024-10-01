# MoranWithSex

Matlab implementations of the codes described in "Fixation probability in a diploid sexually reproducing population", by Zhenyu Shi, Loïc Marrec, and Xiang-Yi Li Richter.

Archived version: TBA

Briefly, we perform stochastic simulations of the Moran model in which we include sexual reproduction and diploidy. 

The codes "Moran_monoecy.m" and "Moran_dioecy.m" simulate the evolutionary dynamics of monoecious and dioecious populations, respectively, and return the probability that a mutation becomes fixed. In order to use the codes, run "pfix = Moran_monoecy(Nrep, delta, h, NAA0, NAB0, NBB0)" or "pfix = Moran_dioecy(Nrep, delta, hF, hM, NAAF0, NABF0, NBBF0, NAAM0, NABM0, NBBM0)" in the command window under Matlab after specifying the parameter values.  

The source code is freely available under the GNU GPLv3 license.

If you find this code useful for your research, please cite the associated reference, "Fixation probability in a diploid sexually reproducing population", by Zhenyu Shi, Loïc Marrec, and Xiang-Yi Li Richter.
