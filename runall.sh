cd COSAMP/
make
cd ../OMPR
make
cd ../OMP
make
cd ../GradeS/
make
cd ..
for i in ../data/*/; do  ./COSAMP/cosamp ${i}X ${i}y ${i}theta; done
for i in ../data/*/; do  ./GradeS/grades ${i}X ${i}y ${i}theta; done
for i in ../data/*/; do  ./OMPR/ompr ${i}X ${i}y ${i}theta; done
for i in ../data/*/; do  ./OMP/omp ${i}X ${i}y ${i}theta; done

