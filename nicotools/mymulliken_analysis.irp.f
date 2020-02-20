program mymulliken_analysis
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  print *, 'Hello world'
  read_wf = .True.
  touch read_wf
  call routine
end

subroutine routine
 implicit none
 call myprint_mulliken
end
