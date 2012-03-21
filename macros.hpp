#include <iostream>
#define PP(x)                  std :: cout << #x << ":" << (x) << std :: endl
#define PP2(x,y)               std :: cout << #x << ',' << #y << ":\t" << (x) << " , " << (y) << std :: endl
#define PP3(x,y,z)             std :: cout << #x << ',' << #y << ',' << #z                              << ":\t" << (x) << " , " << (y) << " , " << (z) << std :: endl
#define PP4(x,y,z,w)           std :: cout << #x << ',' << #y << ',' << #z << ',' << #w                 << ":\t" << (x) << " , " << (y) << " , " << (z) << " , " << w << std :: endl
#define PP5(x,y,z,w,v)         std :: cout << #x << ',' << #y << ',' << #z << ',' << #w << ','          << ":\t" << (x) << " , " << (y) << " , " << (z) << " , " << w << ',' << v << std :: endl
#define unless(x) if(!(x))
#define DYINGWORDS(x) for (int klsdjfslkfj = (x) ? 0 : 1; klsdjfslkfj!=0; klsdjfslkfj--, ({ assert (x); }) )
#define VERYCLOSE(a,b) (1e-07 > fabs((a)-(b)))
#define For(it, container) for( typeof((container).begin()) it = (container).begin(); it != (container).end(); ++it)
#define ELAPSED() (double(clock()) / CLOCKS_PER_SEC)
