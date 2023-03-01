static char help[] = "Reads a PETSc matrix and vector from a file and solves the normal equations.\n\n";

#include <petscksp.h>
#include <petscviewerhdf5.h>

int main(int argc,char **args)
{
  Mat            A,y;            /* matrix, RHS (column vector) */
  Vec            x,r,b;          /* approx solution, RHS (vector), residual */
  PetscErrorCode ierr;
  PetscViewer    fd;              
  char           file[PETSC_MAX_PATH_LEN]="";     /* input file name */
  char           A_name[128]="A",y_name[128]="y",x0_name[128]="x0";  /* name of the matrix, RHS and initial guess */
  PetscInt       its,n,m;
  KSP            ksp;
  MPI_Comm       comm;
  PetscReal      norm, norm1, norm2;
  PetscMPIInt    npe,mype;
  PetscLogDouble tsetup,tsetup1,tsetup2,tsolve,tsolve1,tsolve2;
  PetscViewer    viewer;
  PetscReal      out[3];

  /*
     Initialize PETSC code with the input arguments
  */
  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
  comm = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(comm, &mype);CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm, &npe);CHKERRQ(ierr);
  /*
     Determine files from which we read the linear system
     (matrix, right-hand-side and initial guess vector).
  */
  ierr = PetscOptionsGetString(NULL,NULL,"-f",file,PETSC_MAX_PATH_LEN,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetString(NULL,NULL,"-A_name",A_name,128,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetString(NULL,NULL,"-y_name",y_name,128,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetString(NULL,NULL,"-x0_name",x0_name,128,NULL);CHKERRQ(ierr);

  /* ------    LOAD STAGE    ------ */
  PetscPreLoadBegin(PETSC_FALSE,"Load system");
  /*
     Use HDF5 Viewer to load the matrix and vector from .mat file
  */
  ierr = PetscViewerHDF5Open(comm,file,FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = PetscViewerPushFormat(fd,PETSC_VIEWER_HDF5_MAT);CHKERRQ(ierr);
  /*
     Load the matrix.
  */
  ierr = MatCreate(comm,&A);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)A,A_name);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatLoad(A,fd);CHKERRQ(ierr);
  ierr = MatGetLocalSize(A,&m,&n);CHKERRQ(ierr);
  /*
     Load the RHS vector.
  */
  ierr = MatCreate(comm,&y);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)y,y_name);CHKERRQ(ierr);
  ierr = MatSetFromOptions(y);CHKERRQ(ierr);
  ierr = MatLoad(y,fd);CHKERRQ(ierr);
  ierr = VecCreate(comm,&b);CHKERRQ(ierr);
  ierr = VecSetSizes(b,m,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(b);CHKERRQ(ierr);
	ierr = MatGetColumnVector(y, b, 0);CHKERRQ(ierr);
  /*
     Create initial guess vector (all zeros).
  */
  ierr = VecCreate(comm,&x);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)x,x0_name);CHKERRQ(ierr);
  ierr = VecSetSizes(x,n,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  ierr = VecSet(x,.0);CHKERRQ(ierr);
  /*
     Destroy the viewer for loading the files.
  */
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

  /* ------    SETUP STAGE    ------ */
  PetscPreLoadStage("KspSetUp");
  PetscTime(&tsetup1);
  /* 
     Setup solver 
  */
  ierr = KSPCreate(comm,&ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp, A, A);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  ierr = KSPSetUp(ksp);CHKERRQ(ierr);
  ierr = KSPSetUpOnBlocks(ksp);CHKERRQ(ierr);
  /* 
     Calculate build time
  */
  PetscTime(&tsetup2);
  tsetup = tsetup2 - tsetup1;

  /* ------    SOLVE STAGE    ------ */
  PetscPreLoadStage("KspSolve");
  PetscTime(&tsolve1);
  /* 
     Solve the system 
  */
  ierr = KSPSolve(ksp, b, x);CHKERRQ(ierr);
  /* 
     Calculate solve time
  */
  PetscTime(&tsolve2);
  tsolve = tsolve2 - tsolve1;

  /* ------    CALCULATE NECESSARY OUTPUTS    ------ */
  PetscPreLoadStage("CleanUp");
  ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
  /* 
     Calculate residual norm
  */
  ierr = PetscObjectSetName((PetscObject)x,"x");CHKERRQ(ierr);
  ierr = VecDuplicate(b,&r);CHKERRQ(ierr);
  ierr = MatMult(A,x,r);CHKERRQ(ierr);
  ierr = VecAXPY(r,-1.0,b);CHKERRQ(ierr);
  ierr = VecNorm(r,NORM_2,&norm1);CHKERRQ(ierr);
  /* 
     Normalize residual norm with the norm of RHS
  */
  ierr = VecNorm(b,NORM_2,&norm2);CHKERRQ(ierr);
  norm = norm1/norm2;
  /* 
     Print and write necessary outputs to the file
  */
  ierr = PetscPrintf(comm,"Number of iterations = %d\n",its);CHKERRQ(ierr);
  ierr = PetscPrintf(comm,"Residual norm = %g\n",(double)norm);CHKERRQ(ierr);
  ierr = PetscPrintf(comm,"Build time = %g\n",(double)tsetup);CHKERRQ(ierr);
  ierr = PetscPrintf(comm,"Solve time = %g\n",(double)tsolve);CHKERRQ(ierr);
	ierr = PetscViewerASCIIOpen(comm, "petscOut.txt", &viewer);CHKERRQ(ierr);
  out[0] = (double)tsetup; out[1] = (double)tsolve; out[2] = (double)norm;
	for (PetscInt i = 0; i < 4; ++i) {
		if (i == 3){
			ierr = PetscViewerASCIIPrintf(viewer, "%g", out[i-1]);CHKERRQ(ierr);
		}
		else if (i==2){
			ierr = PetscViewerASCIIPrintf(viewer, "%d ", its);CHKERRQ(ierr);
		}
		else{
			ierr = PetscViewerASCIIPrintf(viewer, "%g ", out[i]);CHKERRQ(ierr);
		}
	}
  ierr = PetscViewerDestroy(&viewer);
  /*
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
  */
  ierr = MatDestroy(&A);CHKERRQ(ierr);   ierr = VecDestroy(&b);CHKERRQ(ierr);
  ierr = VecDestroy(&r);CHKERRQ(ierr);   ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  PetscPreLoadEnd();

  ierr = PetscFinalize();

  return ierr;
}
