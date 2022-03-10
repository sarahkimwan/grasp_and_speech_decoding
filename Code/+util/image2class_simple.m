function [ class ] = image2class_simple( imagename )

%IMAGE2CLASS_SIMPLE returns the name of a grasp if provided a number, or the associated number of the
%grasp if provided a name. 

if ~isnumeric(imagename) 

    if mean(ismember('ExtensionType', imagename)) ==1 
        class =1;

     elseif mean(ismember('InferiorPincer', imagename))==1
        class =2;

    elseif mean(ismember('Lateral', imagename))==1
        class =3;

    elseif mean(ismember('PowerSphere', imagename))==1
        class =4;

    elseif mean(ismember('WritingTripod', imagename))==1
        class =5;

    elseif mean(ismember('MediumWrap', imagename))==1
        class =6;

    elseif mean(ismember('PalmarPinch', imagename))==1
        class =7;

    elseif mean(ismember('Sphere3Finger', imagename))==1 
        class =8;

    elseif mean(ismember('largeDiameter', imagename))==1 
        class =11;

    elseif mean(ismember('Blue', imagename))==1 
        class =12;

    elseif mean(ismember('Red', imagename))==1 
        class =13;
    
    elseif mean(ismember('Green', imagename))==1 
        class =14;

    elseif mean(ismember('Yellow', imagename))==1 
        class =15;

    elseif mean(ismember('Gray', imagename))==1 
        class =16;

    elseif mean(ismember('Brown', imagename))==1 
        class =17;

    elseif mean(ismember('Lateral No Go', imagename))==1
        class =18;

    elseif mean(ismember('WritingTripod No Go', imagename))==1
        class =20;

    elseif mean(ismember('MediumWrap No Go', imagename))==1
        class =21;

    elseif mean(ismember('PalmarPinch No Go', imagename))==1
        class =22;

    elseif mean(ismember('Sphere3Finger No Go', imagename))==1 
        class =23;
            
    elseif mean(ismember('Lateral Image', imagename))==1 
        class =24;

    elseif mean(ismember('WritingTripod Image', imagename))==1 
        class =25;

    elseif mean(ismember('MediumWrap Image', imagename))==1 
        class =26;

    elseif mean(ismember('PalmarPinch Image', imagename))==1 
        class =27;

    elseif mean(ismember('Sphere3Finger Image', imagename))==1 
        class =28;
    
    else
        error([ imagename{1} 'Unknown grasp, add it to list']);
    end
    
elseif(isnumeric(imagename))
    
    if imagename ==1 
     class = 'ExtensionType';

    elseif imagename == 2
     class ='InferiorPincer';

    elseif  imagename == 3
     class = 'Lateral';

    elseif imagename ==4
      class ='PowerSphere'; 

    elseif imagename== 5
       class = 'WritingTripod';

    elseif imagename ==6
       class ='MediumWrap';

    elseif imagename == 7
       class ='PalmarPinch';

    elseif imagename == 8
        class ='Sphere3Finger';

    elseif imagename ==11 
        class ='largeDiameter';

    elseif imagename ==12
        class ='Blue';

    elseif imagename ==13 
        class ='Red';

    elseif imagename ==14
        class ='Green';

    elseif imagename ==15
        class ='Yellow';

    elseif imagename ==16 
        class ='Gray';

    elseif imagename ==17 
        class ='Brown';            

    elseif imagename ==18
        class ='Lateral No Go';

    elseif imagename ==20
        class ='WritingTripod No Go';

    elseif imagename ==21
        class ='MediumWrap No Go';

    elseif imagename ==22
        class ='PalmarPinch No Go';

    elseif imagename ==23
        class ='Sphere3Finger No Go';

    elseif imagename ==24 
        class ='Lateral Image';            

    elseif imagename ==25
        class ='WritingTripod Image';

    elseif imagename ==26
        class ='MediumWrap Image';

    elseif imagename ==27
        class ='PalmarPinch Image';

    elseif imagename ==28
        class ='Sphere3Finger Image';                          
    else
    error([ imagename 'Unknown grasp, add it to list']);
    
    end
    
end 





end

