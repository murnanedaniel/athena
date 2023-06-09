package atlantis.interactions;


import java.awt.*;
import java.awt.geom.*;


/**
 * A panel which selects a parallogram-shaped region in which two
 * sides are parallel to the y-axis and the other two are skewed with
 * respect to the x-axis.
 *
 * @author Charles Loomis
 **/
public class AYSkewSelection extends ASelection {

  final private static int STARTING_WIDTH=25;

  public AYSkewSelection() {
    super(6);
  }

  /**
   * Initialize the control points based in the starting point
   */
  public int init(Point2D.Double p, int key) {
    isValid=false;

    setCenter(hr[0], p.x, p.y-STARTING_WIDTH);
    setCenter(hr[1], p.x, p.y-STARTING_WIDTH);
    setCenter(hr[2], p.x, p.y+STARTING_WIDTH);
    setCenter(hr[3], p.x, p.y+STARTING_WIDTH);
    setCenter(hr[4], p.x, p.y);
    setCenter(hr[5], p.x, p.y);

    region=5;
    return region;
  }

  /**
   * Move the active control point to the point (x,y).
   */
  public void drag(Point2D.Double p, int region, int key) {
    isValid=true;

    double width;

    // Change what is done depending on which control point is active.
    switch(region) {
    case 0:
      width=p.y-hr[4].getCenterY();
      setCenterY(hr[0], hr[4].getCenterY()+width);
      setCenterY(hr[1], hr[5].getCenterY()+width);
      setCenterY(hr[2], hr[5].getCenterY()-width);
      setCenterY(hr[3], hr[4].getCenterY()-width);
      break;

    case 1:
      width=p.y-hr[5].getCenterY();
      setCenterY(hr[0], hr[4].getCenterY()+width);
      setCenterY(hr[1], hr[5].getCenterY()+width);
      setCenterY(hr[2], hr[5].getCenterY()-width);
      setCenterY(hr[3], hr[4].getCenterY()-width);
      break;

    case 2:
      width=p.y-hr[5].getCenterY();
      setCenterY(hr[0], hr[4].getCenterY()-width);
      setCenterY(hr[1], hr[5].getCenterY()-width);
      setCenterY(hr[2], hr[5].getCenterY()+width);
      setCenterY(hr[3], hr[4].getCenterY()+width);
      break;

    case 3:
      width=p.y-hr[4].getCenterY();
      setCenterY(hr[0], hr[4].getCenterY()-width);
      setCenterY(hr[1], hr[5].getCenterY()-width);
      setCenterY(hr[2], hr[5].getCenterY()+width);
      setCenterY(hr[3], hr[4].getCenterY()+width);
      break;

    case 4:
      width=hr[4].getCenterY()-hr[0].getCenterY();
      setCenter(hr[region], p.x, p.y);
      setCenter(hr[0], p.x, p.y-width);
      setCenter(hr[3], p.x, p.y+width);
      break;

    case 5:
      width=hr[4].getCenterY()-hr[0].getCenterY();
      setCenter(hr[region], p.x, p.y);
      setCenter(hr[1], p.x, p.y-width);
      setCenter(hr[2], p.x, p.y+width);
      break;
    }
  }

  public void paint(Graphics2D g) {
    paintStandard(g);
  }

  /**
   * Make the affine transform which corresponds to this skewed region.
   */
  public Point2D.Double[] getCorners() {
    int first=getUpperLeftRegion();

    // Calculate which is the opposite corner.
    int third=(first+2)%4;

    // Now use the cross-product to determine which of the
    // remaining points is the one which keep the path going clockwise.
    int second=(first+1)%4;
    double dx0=hr[third].getCenterX()-hr[first].getCenterX();
    double dy0=hr[third].getCenterY()-hr[first].getCenterY();
    double dx1=hr[second].getCenterX()-hr[first].getCenterX();
    double dy1=hr[second].getCenterY()-hr[first].getCenterY();

    if(dx0*dy1-dy0*dx1>0) second=(first+3)%4;

    return convert(first, second, third);
  }

}
