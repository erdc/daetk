#include "MeshMHFE.h"


MeshMHFE::MeshMHFE(ParameterDatabase& pd):
  nCellsX(pd.i("nxNodes")),
  nCellsY(pd.i("nyNodes")),
  nCellsZ(pd.i("nzNodes")),
  xLeft(pd.r("xLeft")),
  yFront(pd.r("yFront")),
  zBottom(pd.r("zBottom")),
  xRight(pd.r("xRight")),
  yBack(pd.r("yBack")),
  zTop(pd.r("zTop")),
  dX(1.0),dY(1.0),dZ(1.0),
  Xc(0.0),Yc(0.0),Zc(0.0),
  Xf(0.0),Yf(0.0),Zf(0.0)

{
  assert(nCellsX >= 1);
  assert(nCellsY >= 1);
  assert(nCellsZ >= 1);

  dX = (xRight-xLeft)/nCellsX;
  dY = (yBack-yFront)/nCellsY;
  dZ = (zTop-zBottom)/nCellsZ;

}
