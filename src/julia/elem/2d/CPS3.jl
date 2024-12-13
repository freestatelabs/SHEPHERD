""" SHEPHERD Element Library
    CPS3: 3-node linear plane-stress element 
"""

mutable struct CPS3 <: Element 

    num::Integer 
    material::Material 
    thickness::Number

    # Calculated when instantiated 
    centroid::Number
end


function stiffness_CPS3(x1, y1, x2, y2, x3, y3, t, D)

    x12 = x1 - x2 
    x13 = x1 - x3 
    x23 = x2 - x3 
    y12 = y1 - y2 
    y13 = y1 - y3
    y23 = y2 - y3

    detJ = x13*y23 - y13*x23
    Ae = abs(detJ/2)
    B = (1/detJ) * [y23 0 -y13 0 y12 0; 0 -x23 0 x13 0 -x12; -x23 y23 x13 -y13 -x12 y12]

    Ke = t * Ae * transpose(B) * D * B
    return Ke, B
end


