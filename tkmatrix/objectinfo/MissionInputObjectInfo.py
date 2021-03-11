from tirma.objectinfo.ObjectInfo import ObjectInfo

class MissionInputObjectInfo(ObjectInfo):
    """
    Implementation of ObjectInfo to be used to characterize objects which are to be loaded from a csv file and knowing
    the source mission id of the input light curve.
    """
    def __init__(self, mission_id, input_file, initial_mask=None, initial_detrend_period=None, star_info=None,
                 aperture_file=None):
        """
        @param mission_id: the mission identifier. TIC ##### for TESS, KIC ##### for Kepler and EPIC ##### for K2.
        @param input_file: the file to be used for loading the light curve
        @param initial_mask: an array of time ranges provided to mask them into the initial object light curve.
        @param initial_detrend_period: integer value specifying a fixed value for an initial period to be detrended
        @param aperture_file: the file containing 1s and 0s specifying the user selected aperture
        from the initial light curve before processing.
        @param star_info: input star information
        """
        super().__init__(initial_mask, initial_detrend_period, star_info, aperture_file)
        self.id = mission_id
        self.input_file = input_file

    def sherlock_id(self):
        return self.id.replace(" ", "") + "_INP"

    def mission_id(self):
        return self.id
