% matlab batch file mbftsrun
  clear

  fprintf ( 1, '\n' );
  fprintf ( 1, 'tsrun\n' );
  fprintf ( 1, 'Run tsrun locally.\n' );
%
%  BATCH defines the job and sends it for execution.
%
  my_job = batch ( 'tsrun', 'Profile', 'local');
%
%  WAIT pauses the MATLAB session til the job completes.
%
  wait ( my_job );
%
%  DIARY displays any messages printed during execution.
%
  diary ( my_job );
%
%  LOAD makes the script's workspace available.
%
% total = total number of primes.
%
  load ( my_job );

%
%  These commands clean up data about the job we no longer need.
%
  delete ( my_job ); %Use delete() for R2012a or later

  fprintf ( 1, '\n' );
  fprintf ( 1, 'tsrun\n' );
  fprintf ( 1, '  Normal end of execution.\n' );