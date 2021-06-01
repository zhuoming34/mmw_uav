% Testing Java API for controlling mouse pointer
% http://www.java2s.com/Code/JavaAPI/java.awt.event/InputEvent.htm
import java.awt.Robot;
import java.awt.event.*;

mouse = Robot;
screenSize = get(0, 'screensize');
%mouse.mouseMove(0, 0);
%get(0,'PointerLocation')

pause(1)
mouse.mouseMove(555, screenSize(4)-568); % (horizontal,verital)
pause(1)
mouse.mousePress(InputEvent.BUTTON1_MASK);
mouse.mouseRelease(InputEvent.BUTTON1_MASK);
pause(2)
mouse.mouseMove(1087, screenSize(4)-568); 
pause(1)
mouse.mousePress(InputEvent.BUTTON1_MASK);
mouse.mouseRelease(InputEvent.BUTTON1_MASK);