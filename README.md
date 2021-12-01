
*Overview*

Swarm algorithm, PSO (Particle Swarm Optimization) for continous functions.

The algorithm starts with a random set of positions, and then through multiple generations or passes it will converge or settle in on a function.

```console
make libpso.a
make test-mfunc
make pyramid.asc
```

You can alter the test-mfunc.c to use a different math function, such as sinc, bumps, etc.

*Installing meshlab*

```console
sudo apt-get install meshlab
```

*Filtering points and exporting to STL*

Open meshlab

File menu, Import mesh, pyramid.asc

Continuing in MeshLab, use the Filters menu.

There are two areas under that menu that we will be working with.

Remeshing, Simplification, and Reconstruction

And also Normals, Curvatures, and Orientation.

Here are the steps:

1. Simplification, Clustering Decimation. Using 0.1 is a good value for world unit.

2. Compute normals and point sets. (Optionally you may need to flip the normals)

3. Surface Reconstruction: Screened Poisson.

That's it! You can now Export Mesh into STL format for 3d printing, or choose other formats to further work with the completed mesh.

*Total error*

A csv file is generated named pso_generations.csv which has the error value for each generation that is calculated. You can plot this in a spreadsheet program, Python, or with the included R script using R Studio for example.

*Developer notes*

The main process also calculates the swarm particles, so if you specify 8 threads then there is actually the main process plus 7 actual threads (or workers). After each generation, one of the actual threads will calculate the total error. As that is taking place, the other workers will come out of idle if they are in that state and advance to calculating the next generation.

In the process_pso, we wait for all the batch process packs to be filled by the threads.
Then we set the threads to wait for a calculation signal, and one of them does it. We check that the total error value has been updated. Finally we advance to the next generation and signal the threads to continue working.

For additional debugging statements, you can run a small number of particles and generations. Turn on debugging by editing pso.c and adjust the dprintf define from emptyprint to debugprint

![Image of output](https://phrasep.com/~lvecsey/software/libpso/screenshot_bumps2-0_scaled.png)
![Image of output](https://phrasep.com/~lvecsey/software/libpso/screenshot_bumps2-1_scaled.png)
![Image of output](https://phrasep.com/~lvecsey/software/libpso/screenshot_bumps2-2_scaled.png)
![Image of output](https://phrasep.com/~lvecsey/software/libpso/screenshot_bumps2-3.png)
