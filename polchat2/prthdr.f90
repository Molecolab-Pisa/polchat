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
