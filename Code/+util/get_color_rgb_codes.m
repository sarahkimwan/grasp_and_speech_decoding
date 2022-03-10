function [colors] = get_color_rgb_codes(imagename)
%GET_COLOR_RGB_CODES outputs color code for input labels. 

colors = cell(1, length(imagename));
for i =1:length(imagename)
    
     if strcmp('Lateral', imagename{i}) ==1 
        colors{i} = [255, 176, 0]/255;
    elseif strcmp('WritingTripod',  imagename{i}) ==1 
        colors{i} = [100, 143, 255]/255;

    elseif strcmp('MediumWrap',  imagename{i}) ==1 
        colors{i} = [254, 97, 0]/255;

    elseif strcmp('PalmarPinch',  imagename{i}) ==1 
        colors{i} = [120, 94, 240]/255;

    elseif strcmp('Sphere3Finger',  imagename{i}) ==1 
        colors{i} = [220, 38, 127]/255;

    elseif strcmp('Lateral No Go',  imagename{i}) ==1 
        colors{i} = util.rgb('Tomato');
        
    elseif strcmp('WritingTripod No Go',  imagename{i}) ==1 
        colors{i} = util.rgb('LightGreen');
        
    elseif strcmp('MediumWrap No Go',  imagename{i}) ==1 
        colors{i} = util.rgb('LightCyan');
        
    elseif strcmp('PalmarPinch No Go',  imagename{i}) ==1 
        colors{i} = util.rgb('Yellow');
        
    elseif strcmp('Sphere3Finger No Go',  imagename{i}) ==1 
        colors{i} = util.rgb('Violet');
 
     elseif strcmp('Lateral Image', imagename{i}) ==1 
        colors{i} = util.rgb('Green');
        
    elseif strcmp('WritingTripod Image',  imagename{i}) ==1 
        colors{i} = util.rgb('LawnGreen');
        
    elseif strcmp('MediumWrap Image',  imagename{i}) ==1 
        colors{i} = util.rgb('MediumAquamarine');
        
    elseif strcmp('PalmarPinch Image',  imagename{i}) ==1 
        colors{i} = util.rgb('SeaGreen');
        
    elseif strcmp('Sphere3Finger Image',  imagename{i}) ==1 
        colors{i} = util.rgb('ForestGreen'); 
                
        
    elseif  strcmp('Action',  imagename{i}) ==1 
        colors{i} = util.rgb('SteelBlue');       
    elseif strcmp('Cue',  imagename{i}) ==1 
        colors{i} = util.rgb('Orange');  
    elseif strcmp('Image Cue',  imagename{i}) ==1 
        colors{i} = util.rgb('Orange');  
    elseif strcmp('CuePhase',  imagename{i}) ==1 
        colors{i} = util.rgb('Orange');  
    elseif strcmp('Both',  imagename{i}) ==1 
        colors{i} = util.rgb('LightSlateGray');       
    elseif strcmp('PMV',  imagename{i}) ==1 
        colors{i} = [0 1 0];
        
    elseif strcmp('SMG',  imagename{i}) ==1 
        colors{i} = [0 0 1]; 
    elseif strcmp('S1X',  imagename{i}) ==1 
        colors{i} = util.rgb('Magenta'); 
    
    elseif  strcmp('ITI',  imagename{i}) ==1 
        colors{i} = util.rgb('Silver');    
        
    elseif  strcmp('Delay',  imagename{i}) ==1 
        colors{i} = util.rgb('Cyan'); 
        
    elseif  strcmp('MotorImagery',  imagename{i}) ==1 
        colors{i} = util.rgb('Green'); 
    
    elseif  strcmp('Speaking',  imagename{i}) 
        colors{i} = util.rgb('Blue'); 
        
    elseif  strcmp('Go Trial',  imagename{i}) 
        colors{i} = util.rgb('Green'); 
        
    elseif  strcmp('No Go Trial',  imagename{i}) 
        colors{i} = util.rgb('Red'); 
        
    elseif  strcmp('Blue',  imagename{i}) 
        colors{i} = util.rgb('Blue'); 
        
    elseif  strcmp('Green',  imagename{i}) 
        colors{i} = util.rgb('Green'); 
        
    elseif  strcmp('Yellow',  imagename{i}) 
        colors{i} = util.rgb('Fire'); 
      
    elseif  strcmp('Brown',  imagename{i}) 
        colors{i} = util.rgb('Brown'); 
          
     elseif  strcmp('Gray',  imagename{i}) 
        colors{i} = util.rgb('SlateGray'); 
         
     elseif  strcmp('Grasps',  imagename{i}) 
        colors{i} = util.rgb('Blue'); 
        
     elseif  strcmp('Colors',  imagename{i}) 
        colors{i} = util.rgb('Red'); 
        
     elseif  strcmp('Colors',  imagename{i}) 
        colors{i} = util.rgb('Red'); 
    
     elseif strcmp('ImageCue', imagename{i})
         colors{i} = util.rgb('Green');
         
      elseif strcmp('Image', imagename{i})
         colors{i} = util.rgb('Green');
         
     elseif strcmp('Training', imagename{i})
     colors{i} = util.rgb('Green');
     
     elseif strcmp('Testing', imagename{i})
     colors{i} = util.rgb('Blue');
     
     elseif strcmp('Motor Imagery', imagename{i}) 
        colors{i} = [0 0 1];      
        
    elseif strcmp('Spoken Grasps', imagename{i})
        colors{i} = [144 144 144]/255;
        
       
     elseif strcmp('Spoken Colors', imagename{i})
        colors{i} = [41 204 230]/255;
        
    else  
        error(['Condition ' imagename{i} ' not present, add color'])
        
    end 
end 

end

