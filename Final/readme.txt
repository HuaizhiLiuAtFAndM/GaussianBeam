polarity.py
    Class vector
        This class defines vectors with 2 parameters: x value and y value (Both value are complex)
    Based on vector class, several common unit vactor used to represent light polarization is defined:
        e_x
        e_y
        e_L
        e_R
        etc.
    
BeamProject.py
    Class Screen
        This class defines a screem of a given z value (distance to source)
        Value of beam is represented by a 2D array of vectors defined in polarity.py
    Class Beam
        This class defines a light beam represented by its basic_parameter (for example, p and l for Laguerre Gaussian Beam), and a expression
        The basicParameter is an 1D array.
        The expression must be a callable with 6 parameters in following order.
            basicParameter, x, y, z, Rho, Theta
            where z is distance to source corresponding to the Screen
            x and y are from cartisan coordinate
            Rho and Theta are from cylindrical coordinate
    Func Render
        This function calculate the beam expression and plot it to the sceen
    Func Export 
        This function get values form the screen and display it or save as a png file.
        polarization mode might be represent as linear or elliptical 

LaguerreGaussianBeam.py
    The file defines a sub class of Beam.
    example shown is radial polarization and azimuthal polarization

