import scipy.optimize as optimize

file = open("out","r")
f_act = file.readline()
f_target = 1244.51

f = lambda f_act,f_target : f_act - f_target

