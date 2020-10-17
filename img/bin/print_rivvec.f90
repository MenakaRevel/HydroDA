      program print_fd
! ===============================================
      implicit none
! input
      character*128      ::  vectxt
      integer            ::  type, level
      character*128      ::  buf
! vars
      real               ::  lon1, lat1, lon2, lat2
      real               ::  uparea
      integer            ::  level_this
! ===============================================
      call getarg(1,vectxt)
      call getarg(2,buf)
      read(buf,*) type
      call getarg(3,buf)
      read(buf,*) level

      open(10,file=vectxt,form='formatted')
 1000 continue
        read(10,*,end=9999) lon1, lat1, lon2, lat2, uparea

        if( type==1 )then
          level_this=int( log10(uparea/100.)*2. )
          level_this=max(level_this,1)
          level_this=min(level_this,10)

          if( level_this==level .and. lon2>-999 ) then
            write(6,*) lon1, lat1, uparea, lon2, lat2
          endif

        elseif( type==2 )then
          write(6,*) lon1, lat1

        endif

      goto 1000
 9999 continue

      close(10)

      end program print_fd
