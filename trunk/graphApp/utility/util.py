"""Miscellaneous functions used here and there"""

def numToText(num, precision, mode):
    """Given a number and the precision, return it as text. Makes 
    sure the correct number of places after the decimal point are shown.
    
    """

    # precision refers to the number of places after to decimal place to show.
    # if a number has too few numbers, it is padded with zeros.
    # if a number has more numbers after the decimal place than the 
    # precision specified, then the argument 'mode' determines what happens.
    # mode can be 'round', 'trunc', or 'extend'
    # mode = 'round' means to round the number
    # mode = 'trunc' means to truncate the number
    # mode = 'extend' means to show all digits of the number, even if it
    #         is more than the precision specified

    # if mode is not understood, default to round
    if mode not in ['round', 'trunc', 'extend']:
        mode = 'round'

    s = str(float(num))
    left, dot, right = s.partition('.')
    if precision == 0:
        # don't show any decimal places. 
        # see if it's already a whole number
        if dot == '' or right == '':
            return left
        # action to do is determined by mode
        if mode == 'round':
            return str(int(round(float(num), 0)))
        if mode == 'trunc':
            return left
        if mode == 'extend':
            # strip out any unnecessary zeros
            r = ''
            hit = False
            for i in range(len(right)-1, -1, -1):
                if not hit and right[i] != '0':
                    hit = True
                if hit:
                    r = right[i] + r
            if r == '':
                return left
            else:
                return left + dot + r
    else:
        # need to show some decimal places
        length = len(right)
        if precision == length:
            return left + dot + right
        elif precision > length:
            # pad it with zeros
            right = right + '0'*(precision-length)
            return left + dot + right
        else:
            # there are more digits than precision
            if mode == 'round':
                return str(round(float(num), precision))
            if mode == 'trunc':
                return left + dot + right[:precision]
            if mode == 'extend':
                return left + dot + right
