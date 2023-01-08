function iCorner = find_corner(x,y,whichCorner)

%FINDCORNER find the corner of a rectangle-like shape.
% ICORNER = FINDCORNER(X,Y,'corner') returns index of X and Y at which the
% corner of a rectangle-like shape is located. The corner is specified by
% the string 'corner' (= 'lowerleft', 'upperright', etc.). 

% J.W. Lankhaar, 20 September 2006.
% Revisions:
% 09-Jan-07 Distance relative to bounding box.

% Find center of the rectangle (center of outer bounding box).
xCenter = (max(x)+min(x))/2;
yCenter = (max(y)+min(y))/2;

% Find x and y width of bounding box.
xWidth = max(x)-min(x);
yWidth = max(y)-min(y);

% Find distance to center (relative to bounding box dimensions).
d = sqrt(((x-xCenter)/xWidth).^2+((y-yCenter)/yWidth).^2);

% Find the appropriate segment of the shape.
switch lower(whichCorner)
	case 'upperleft'
		iSegment = find(x < xCenter & y > yCenter);
	case 'upperright'
		iSegment = find(x > xCenter & y > yCenter);
	case 'lowerleft'
		iSegment = find(x < xCenter & y < yCenter);
	case 'lowerright'
		iSegment = find(x > xCenter & y < yCenter);
end

% Find the corner (maximum distance to center).
[xx,iMax] = max(d(iSegment));
iCorner = iSegment(iMax);