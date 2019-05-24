obj/LeeSuzuki.o : src/LeeSuzuki.F90 obj/Profiler.o obj/LinAlgLib.o 
obj/UMOA.o : src/UMOA.F90 obj/LeeSuzuki.o obj/Iteration.o obj/Operators.o obj/LinAlgLib.o 
obj/LinAlgLib.o : submodule/LinAlgf90/src/LinAlgLib.f90 obj/MatVecComplex.o obj/MatVecDouble.o obj/MatVecSingle.o obj/MatrixComplex.o obj/MatrixDouble.o obj/MatrixSingle.o obj/VectorComplex.o obj/VectorDouble.o obj/VectorSingle.o obj/SingleDoubleComplex.o obj/LinAlgParameters.o 
obj/MatVecDouble.o : submodule/LinAlgf90/src/MatVecDouble.f90 obj/MatrixDouble.o obj/VectorDouble.o obj/LinAlgParameters.o 
obj/VectorSingle.o : submodule/LinAlgf90/src/VectorSingle.f90 obj/LinAlgParameters.o 
obj/LinAlgParameters.o : submodule/LinAlgf90/src/LinAlgParameters.f90 
obj/VectorComplex.o : submodule/LinAlgf90/src/VectorComplex.f90 obj/LinAlgParameters.o 
obj/MatrixSingle.o : submodule/LinAlgf90/src/MatrixSingle.f90 obj/VectorSingle.o obj/LinAlgParameters.o 
obj/MatrixDouble.o : submodule/LinAlgf90/src/MatrixDouble.f90 obj/VectorDouble.o obj/LinAlgParameters.o 
obj/MatVecSingle.o : submodule/LinAlgf90/src/MatVecSingle.f90 obj/MatrixSingle.o obj/VectorSingle.o obj/LinAlgParameters.o 
obj/VectorDouble.o : submodule/LinAlgf90/src/VectorDouble.f90 obj/LinAlgParameters.o 
obj/MatrixComplex.o : submodule/LinAlgf90/src/MatrixComplex.f90 obj/VectorComplex.o obj/LinAlgParameters.o 
obj/SingleDoubleComplex.o : submodule/LinAlgf90/src/SingleDoubleComplex.f90 obj/MatrixComplex.o obj/VectorComplex.o obj/MatrixDouble.o obj/VectorDouble.o obj/MatrixSingle.o obj/VectorSingle.o 
obj/MatVecComplex.o : submodule/LinAlgf90/src/MatVecComplex.f90 obj/MatrixComplex.o obj/VectorComplex.o obj/LinAlgParameters.o 
obj/MyLibrary.o : submodule/HartreeFock/src/MyLibrary.f90 obj/ClassSys.o 
obj/ClassSys.o : submodule/HartreeFock/src/ClassSys.f90 
obj/ThreeBodyInteraction.o : submodule/HartreeFock/src/ThreeBodyInteraction.F90 obj/TwoBodyOperator.o obj/LinAlgLib.o obj/ClassSys.o obj/MyLibrary.o obj/ThreeBodyModelSpace.o obj/SingleParticleState.o obj/Profiler.o 
obj/Operators.o : submodule/HartreeFock/src/Operators.F90 obj/Profiler.o obj/DefineOperators.o obj/ThreeBodyInteraction.o obj/ThreeBodyOperator.o obj/TwoBodyOperator.o obj/OneBodyOperator.o obj/ModelSpace.o 
obj/Profiler.o : submodule/HartreeFock/src/Profiler.F90 obj/MPIFunction.o obj/ClassSys.o 
obj/MBPT.o : submodule/HartreeFock/src/MBPT.F90 obj/MyLibrary.o obj/Profiler.o obj/StoreCouplings.o obj/Operators.o obj/ModelSpace.o 
obj/ThreeBodyModelSpace.o : submodule/HartreeFock/src/ThreeBodyModelSpace.F90 obj/MyLibrary.o obj/LinAlgLib.o obj/SingleParticleState.o 
obj/MPIFunction.o : submodule/HartreeFock/src/MPIFunction.F90 
obj/ModelSpace.o : submodule/HartreeFock/src/ModelSpace.F90 obj/MyLibrary.o obj/Profiler.o obj/ClassSys.o obj/ThreeBodyModelSpace.o obj/TwoBodyModelSpace.o obj/OneBodyModelSpace.o obj/SingleParticleState.o 
obj/StoreCouplings.o : submodule/HartreeFock/src/StoreCouplings.F90 obj/MyLibrary.o obj/Profiler.o 
obj/OneBodyOperator.o : submodule/HartreeFock/src/OneBodyOperator.F90 obj/ClassSys.o obj/DefineOperators.o obj/MyLibrary.o obj/OneBodyModelSpace.o obj/LinAlgLib.o 
obj/SingleParticleState.o : submodule/HartreeFock/src/SingleParticleState.F90 obj/ClassSys.o 
obj/TwoBodyOperator.o : submodule/HartreeFock/src/TwoBodyOperator.F90 obj/Profiler.o obj/ClassSys.o obj/DefineOperators.o obj/MyLibrary.o obj/OneBodyOperator.o obj/TwoBodyModelSpace.o obj/LinAlgLib.o 
obj/OneBodyModelSpace.o : submodule/HartreeFock/src/OneBodyModelSpace.F90 obj/SingleParticleState.o 
obj/DefineOperators.o : submodule/HartreeFock/src/DefineOperators.F90 obj/MyLibrary.o 
obj/TwoBodyModelSpace.o : submodule/HartreeFock/src/TwoBodyModelSpace.F90 obj/MyLibrary.o obj/SingleParticleState.o 
obj/HartreeFock.o : submodule/HartreeFock/src/HartreeFock.F90 obj/MyLibrary.o obj/Profiler.o obj/Operators.o obj/LinAlgLib.o 
obj/ThreeBodyOperator.o : submodule/HartreeFock/src/ThreeBodyOperator.F90 obj/ClassSys.o obj/ThreeBodyInteraction.o obj/MyLibrary.o obj/TwoBodyOperator.o obj/OneBodyOperator.o obj/ThreeBodyModelSpace.o obj/LinAlgLib.o 
obj/Iteration.o : submodule/Iteration/src/Iteration.F90 obj/LinAlgLib.o 
obj/WriteOperator.o : main/WriteOperator.F90 obj/Profiler.o obj/ClassSys.o obj/Operators.o obj/UMOAInput.o 
obj/UMOAInput.o : main/UMOAInput.F90 obj/ClassSys.o 
obj/UMOAMain.o : main/UMOAMain.F90 obj/WriteOperator.o obj/UMOA.o obj/HartreeFock.o obj/Operators.o obj/ModelSpace.o obj/UMOAInput.o obj/Profiler.o 
