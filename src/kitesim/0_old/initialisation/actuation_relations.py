
# actuation relations
def up_to_ld(u_p,delta_ld_used):
    '''
    input:  u_p         = [0-1]   (power-setting)
            delta_ld    = 8 or 13 [%] (depower-tape extension, in percentage)
    output: depower-tape extension [mm]
    '''
    if delta_ld_used == 8:
        depower_tape_max = 1.482
    elif delta_ld_used == 13:
        depower_tape_max = 1.722
    else:
        raise Exception('delta_ld wrong value, should be: 8 or 13')

    depower_tape_0 = 1.098 #minimum depower tape length

    return 0.5*(1-u_p)*(depower_tape_max-depower_tape_0)

### TURNING 
# def us_to_ls(u_s):
#     #u_s = -0.5 -> min_steer and u_s = 0.5 -> max_steer
#     # negative u_s means steering to the right!
#     min_steer, max_steer = -40,40  # Defining the min/max steering numbers observed from the measurements (averaged)
#     ratio = (max_steer-min_steer)*u_s #This will be number between 0 and (-)80
#     ls = delta_ls_max*(ratio/100)
#     ls_left,ls_right = ls_0 - ls, ls_0 + ls
#     return ls_left,ls_right
