function intersection = HeadingIntersection(p1, heading1, p2, heading2)

    ax1 = cos(heading1);
    ay1 = sin(heading1);
    bx1 = p1(1);
    by1 = p1(2);

    ax2 = cos(heading2);
    ay2 = sin(heading2);
    bx2 = p2(1);
    by2 = p2(2);

    t1 = (by2 + ay2 / ax2 * (bx1 - bx2) - by1) / (ay1 - ay2 * ax1 / ax2)

    intersection = [bx1 + ax1 * t1
                    by1 + ay1 * t1];