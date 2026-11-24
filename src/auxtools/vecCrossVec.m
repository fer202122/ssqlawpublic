function output = vecCrossVec(vec1, vec2)
%VECCROSSVEC Cross produc of two vector time-series (i.e., [x y z] * [u v w])
    output = cross(vec1, vec2, 2);
end