To make `pendulum` 

* you should have installed capd dynamical systems library (see `...git/mrs2/companions/README.md`)
* modify the right path to `...git/mrs2/companions/capd-capdDynSys-4.2.153/bin` directory.

```%sh
$ make pendulum
g++ -O2 `.../git/mrs2/companions/capd-capdDynSys-4.2.153/bin/capd-config --cflags` pendulum.cpp `.../git/mrs2/companions/capd-capdDynSys-4.2.153/bin/capd-config --libs` -o pendulum

$ ./pendulum
Syntax: pendulum theta1 theta2 Nt g

$ ./pendulum 1 2 10 9.8
+0.0000000000e+00 -0.0000000000e+00 +0.0000000000e+00 -0.0000000000e+00 +0.0000000000e+00 +1.7453292520e-02 +1.7453292520e-02 +3.3906585040e-02 +3.5906585040e-02 
+1.0000000001e-01 +1.7619410852e-03 +1.7638301917e-03 +3.1369504882e-03 +3.3217889054e-03 +1.7961760201e-02 +1.8016412016e-02 +2.6347744526e-02 +2.8326159811e-02 
+2.0000000001e-01 +3.6165519172e-03 +3.6301394323e-03 +4.9400705912e-03 +5.2302256538e-03 +1.9174028121e-02 +1.9364266138e-02 +7.5896542181e-03 +9.3429830246e-03 
+3.0000000001e-01 +5.5865851678e-03 +5.6266124014e-03 +4.6769732201e-03 +4.9493919966e-03 +2.0069382795e-02 +2.0411945687e-02 -1.4605503914e-02 -1.2680056260e-02 
+4.0000000001e-01 +7.5807674299e-03 +7.6609629554e-03 +2.4143110874e-03 +2.7125458125e-03 +1.9459342758e-02 +1.9955819556e-02 -3.0510943832e-02 -2.8770712481e-02 
+5.0000000001e-01 +9.4038898059e-03 +9.5367220653e-03 -7.9779842545e-04 -5.6681908722e-04 +1.6548934493e-02 +1.7185441968e-02 -3.3817521309e-02 -3.1884065544e-02 
+6.0000000001e-01 +1.0835004235e-02 +1.0995544447e-02 -3.6407572653e-03 -3.4255986273e-03 +1.1459626906e-02 +1.1909978060e-02 -2.3022984559e-02 -2.1101767395e-02 
+7.0000000001e-01 +1.1653274929e-02 +1.1856382346e-02 -4.9260735889e-03 -4.6355838901e-03 +4.6351267054e-03 +5.2609293908e-03 -2.7922363095e-03 -1.1906545216e-03 
+8.0000000001e-01 +1.1807652354e-02 +1.1986560912e-02 -4.0517198480e-03 -3.7648424910e-03 -2.0955598249e-03 -1.9464256000e-03 +1.8223697878e-02 +1.9311839079e-02 
+9.0000000001e-01 +1.1291967008e-02 +1.1484952421e-02 -1.4241745362e-03 -1.1798027609e-03 -8.0514203917e-03 -7.6928521993e-03 +3.0496945103e-02 +3.2321084038e-02 
+1.0000000000e+00 +1.0327384094e-02 +1.0452099140e-02 +1.8690020185e-03 +1.9703575588e-03 -1.1912989530e-02 -1.1565141227e-02 +2.9703795498e-02 +3.1567915882e-02 
```