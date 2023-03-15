Simulations for testing the stability:

ANAL tests:
Testing the analytical nr input to see how big of an expansion is needed to 
properly simulate with cmb:
[25 50 50 5] - 0 ok
[5 50 50 25] - 1 ok
[25 50 50 25] - 2 ok
[5 100 100 5] - 3 ok
[25 100 100 5] - 4 ok
[50 100 100 5] - 5 ok
[50 100 100 25] - 6 ok

BDR tests:
Testing how big should the "stretching" zone above and below the CMB should be:
[3480, 3550] - 0 error (too tight)
[3479, 3551] - 1 error (too tight)
[3445, 3585] - 2 ok
[3445, 3620] - 3 ok
[3410, 3585] - 4 ok
[3410, 3620] - 5 fail (I basically copied BDR0)

ANA tests:
Testing the control depths for analytical nr:
[0., 2821e3, 2891e3, 6371e3] - 0 ok
[0., 2820e3, 2892e3, 6371e3] - 1 ok
[0., 2810e3, 2900e3, 6371e3] - 2 ok
[0., 2750e3, 2900e3, 6371e3] - 3 ok
[0., 2750e3, 2950e3, 6371e3] - 4 ok

CST tests:
Testing with constant nr everywhere.
5 - 0 ok
10 - 1 ok
25 - 2 ok
50 - 3 ok
100 - 4 ok