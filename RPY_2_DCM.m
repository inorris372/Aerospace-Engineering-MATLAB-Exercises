function dcm = RPY_2_DCM(roll, pitch, yaw)
  cR = cos(roll);  sR = sin(roll);
  cP = cos(pitch); sP = sin(pitch);
  cY = cos(yaw);   sY = sin(yaw);

  dcm = [cP*cY          cP*sY          -sP;
         sR*sP*cY-cR*sY sR*sP*sY+cR*cY sR*cP;
         cR*sP*cY+sR*sY cR*sP*sY-sR*cY cR*cP];
end

