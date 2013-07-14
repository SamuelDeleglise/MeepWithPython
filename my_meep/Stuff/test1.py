class epsilon(CallbackMatrix2D):
    def __init__(self,resolution):
        CallbackMatrix2D.__init__(self)
        master_printf("Creating the material matrix....\n")
        self.meep_resolution = int(resolution)
        eps_matrix = numpy.zeros([10.0 * resolution, 10.0 * resolution], dtype = float)
        len_x = eps_matrix.shape[0]
        len_y = eps_matrix.shape[1]
        for x in range(0, len_x):
            for y in range(0, len_y):
                if ((y >= 4.0) and (y <= 5.0)):
                     eps_matrix[x,y] = 12.0;
                else:
                     eps_matrix[x,y] = 1.0;
        master_printf("Setting the material matrix...\n")
        self.set_matrix_2D(eps_matrix)
        self.stored_eps_matrix = eps_matrix #to prevent the garbage collector from cleaning up the matrix...
        master_printf("MeepMaterial object initialized.\n")
