! prthdr.f90:      A Polarisation consistent charge-fitting tool 
!                  A Molecolab Tool www.molecolab.dcci.unipi.it/tools
!
! Copyright (C) 2014, 2015, 2016, 2017
!   S. Caprasecca, C. Curutchet, S. Jurinovich, B. Mennucci
!
! This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
! A copy of the GNU General Public License can be found in LICENSE or at
!   <http://www.gnu.org/licenses/>.
!
subroutine PrtHdr

  use constants

 1000 format(/,&
     '                                                                           ',/,&
     '                                       mm                                  ',/,&
     '                                    mMMm                                   ',/,&
     '                                  mMMMMm         m                         ',/,&
     '                                 mMMMMm          mMm                       ',/,&
     '                                 mMMMMm          mMm                       ',/,&
     '                                 mMMMMMm        mMMm                       ',/,&
     '                                 MMMMMMMMMMMMMMMMMMm                       ',/,&
     '                                mMMMMMMMMMMMMMMMMMm                        ',/,&
     '       __  ___      __    ____________      __MMMm     __                  ',/,&
     '      /  |/  /___  / /   / ____/ ____/___  / /  ____ _/ /_                 ',/,&
     '     / /|_/ / __ \/ /   / __/ / /   / __ \/ /  / __ `/ __ \                ',/,&
     '    / /  / / /_/ / /___/ /___/ /___/ /_/ / /__/ /_/ / /_/ /                ',/,&
     '   /_/  /_/\__________/_____/\____/_____/_____|__,_/_.___/                 ',/,&
     '           /_  __/ __ \/ __ \/ /  / ___/                                   ',/,&
     '            / / / / / / / / / /   \__ \                                    ',/,&
     '           / / / /_/ / /_/ / /___ __/ /                                    ',/,&
     '          /_/  \____/\____/_____/____/                                     ',/,&
     '            mMMMMMMMMMMMMMMMm                                              ',/,&
     '          mMMMMMMMMMMMMMMMm                                                ',/,&
     '        mMMMMMMMMMMMMMMMMM   + ------------------------------------ +      ',/,&
     '       mMMMMMMMMMMMMMMMMm    |            P O L C H A T             |      ',/,&
     '      mMMMMMMMMMMMMMMMMMm    + ------------------------------------ +      ',/,&
     '      mMMMMm       mMMMMMm   | Stefano Caprasecca                   |      ',/,&
     '      mMMMm       mMMMMMMm   | Carles Curutchet                     |      ',/,&
     '       mMm       mMMMMMMm    | Sandro Jurinovich                    |      ',/,&
     '        m       mMMMMMMm     | Benedetta Mennucci         ver ',A5,' |     ',/,&
     '               mMMMMMm       |          www.molecolab.dcci.unipi.it |     ',/,&
     '                             + ------------------------------------ +      ',/,/,&
     '  Please cite this tool as: ',/, &
     '    PolChat: A polarisation-consistent charge fitting tool. ',/,&
     '    S. Caprasecca, C. Curutchet, S. Jurinovich, B. Mennucci. ',/,&
     '    Molecolab Tools. 2014 www.molecolab.dcci.unipi.it/tools',/)

  write(IOut,1000) Version

  return
end subroutine
