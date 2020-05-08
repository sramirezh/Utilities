#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 09:01:03 2020
based on based.py, LAMMPS.py and LAMMPSparser.py from mdanalysis
@author: simon
"""
import numpy as np
import os
import six
import pandas as pd
from io import StringIO

# TODO add the value of the timestep
class timestep(object):
    def __init__(self, frame, n_chunks, data):
        """
        data is a dataframe
        """
        self.frame = frame
        self._n_chunks = n_chunks
        self.data = data
        self.header = []
                
    def __len__(self):
        return self.n_chunks

    def __iter__(self):
        """Iterate over chunks

        ``for chunk in timestep``

            iterate of the chunks
        """
        for i in range(self.n_chunks):
            yield self[i]
    


class FrameIteratorBase(object):
    """
    Base iterable over the frames of a trajectory.

    A frame iterable has a length that can be accessed with the :func:`len`
    function, and can be indexed similarly to a full trajectory. When indexed,
    indices are resolved relative to the iterable and not relative to the
    trajectory.

    .. versionadded:: 0.19.0

    """
    def __init__(self, trajectory):
        self._trajectory = trajectory

    def __len__(self):
        raise NotImplementedError()

    @staticmethod
    def _avoid_bool_list(frames):
        if isinstance(frames, list) and frames and isinstance(frames[0], bool):
            return np.array(frames, dtype=bool)
        return frames

    @property
    def trajectory(self):
        return self._trajectory


class FrameIteratorSliced(FrameIteratorBase):
    """
    Iterable over the frames of a trajectory on the basis of a slice.

    Parameters
    ----------
    trajectory: ProtoReader
        The trajectory over which to iterate.
    frames: slice
        A slice to select the frames of interest.

    See Also
    --------
    FrameIteratorBase

    .. versionadded:: 0.19.0

    """
    def __init__(self, trajectory, frames):
        # It would be easier to store directly a range object, as it would
        # store its parameters in a single place, calculate its length, and
        # take care of most the indexing. Though, doing so is not compatible
        # with python 2 where xrange (or range with six) is only an iterator.
        super(FrameIteratorSliced, self).__init__(trajectory)
        self._start, self._stop, self._step = trajectory.check_slice_indices(
            frames.start, frames.stop, frames.step,
        )

    def __len__(self):
        return range_length(self.start, self.stop, self.step)

    def __iter__(self):
        for i in range(self.start, self.stop, self.step):
            yield self.trajectory[i]
        self.trajectory.rewind()

    def __getitem__(self, frame):
        if isinstance(frame, numbers.Integral):
            length = len(self)
            if not -length < frame < length:
                raise IndexError('Index {} is out of range of the range of length {}.'
                                 .format(frame, length))
            if frame < 0:
                frame = len(self) + frame
            frame = self.start + frame * self.step
            return self.trajectory._read_frame_with_aux(frame)
        elif isinstance(frame, slice):
            step = (frame.step or 1) * self.step
            if frame.start is None:
                if frame.step is None or frame.step > 0:
                    start = self.start
                else:
                    start = self.start + (len(self) - 1) * self.step
            else:
                start = self.start + (frame.start or 0) * self.step
            if frame.stop is None:
                if frame.step is None or frame.step > 0:
                    last = start + (range_length(start, self.stop, step) - 1) * step
                else:
                    last = self.start
                stop = last + np.sign(step)
            else:
                stop = self.start + (frame.stop or 0) * self.step

            new_slice = slice(start, stop, step)
            frame_iterator = FrameIteratorSliced(self.trajectory, new_slice)
            # The __init__ of FrameIteratorSliced does some conversion between
            # the way indices are handled in slices and the way they are
            # handled by range. We need to overwrite this conversion as we
            # already use the logic for range.
            frame_iterator._start = start
            frame_iterator._stop = stop
            frame_iterator._step = step
            return frame_iterator
        else:
            # Indexing with a lists of bools does not behave the same in all
            # version of numpy.
            frame = self._avoid_bool_list(frame)
            frames = np.array(list(range(self.start, self.stop, self.step)))[frame]
            return FrameIteratorIndices(self.trajectory, frames)

    @property
    def start(self):
        return self._start

    @property
    def stop(self):
        return self._stop

    @property
    def step(self):
        return self._step


class FrameIteratorAll(FrameIteratorBase):
    """
    Iterable over all the frames of a trajectory.

    Parameters
    ----------
    trajectory: ProtoReader
        The trajectory over which to iterate.

    See Also
    --------
    FrameIteratorBase

    .. versionadded:: 0.19.0

    """
    def __init__(self, trajectory):
        super(FrameIteratorAll, self).__init__(trajectory)

    def __len__(self):
        return self.trajectory.n_frames

    def __iter__(self):
        return iter(self.trajectory)

    def __getitem__(self, frame):
        return self.trajectory[frame]


class FrameIteratorIndices(FrameIteratorBase):
    """
    Iterable over the frames of a trajectory listed in a sequence of indices.

    Parameters
    ----------
    trajectory: ProtoReader
        The trajectory over which to iterate.
    frames: sequence
        A sequence of indices.

    See Also
    --------
    FrameIteratorBase
    """
    def __init__(self, trajectory, frames):
        super(FrameIteratorIndices, self).__init__(trajectory)
        self._frames = []
        for frame in frames:
            if not isinstance(frame, numbers.Integral):
                raise TypeError("Frames indices must be integers.")
            frame = trajectory._apply_limits(frame)
            self._frames.append(frame)
        self._frames = tuple(self._frames)

    def __len__(self):
        return len(self.frames)

    def __iter__(self):
        for frame in self.frames:
            yield self.trajectory._read_frame_with_aux(frame)

    def __getitem__(self, frame):
        if isinstance(frame, numbers.Integral):
            frame = self.frames[frame]
            return self.trajectory._read_frame_with_aux(frame)
        else:
            frame = self._avoid_bool_list(frame)
            frames = np.array(self.frames)[frame]
            return FrameIteratorIndices(self.trajectory, frames)

    @property
    def frames(self):
        return self._frames





class chunk_reader(object):
    
    def __init__(self, filename):
        self.filename = filename
#        self._read_first_frame()     #An underscore (_) at the beginning is used to denote private variables in Python
#        self._reopen()
        self._get_initial()
        self._get_position_frames()
        self.n_frames = len(self.offsets)-1
        self.frames = []
#        self._read_next_timestep()
#        root, ext = os.path.splitext(self.filename)
#        self._cache = {}
        
        
#    def _reopen(self):
#        self.close()
#        self._file = open(self.filename)
##        self.ts = timestep(self.n_chunks) # Creates an element of the class timestep
#        self.ts.frame = -1
    
    def close(self):
        if hasattr(self, '_file'):
            self._file.close()
    
    def _get_initial(self):
        """
        reads the header and estimates the number of chunks
        """
        header_lines = 0
        f = open( self.filename)
        while(f.read(1)=='#'):
              header_lines += 1
              self.header = f.readline().split()
        self._first_line = f.tell()-1      
        f.seek(0)
        f.seek(self._first_line)
        self.n_chunks = int(f.readline().split()[1])  # Timestep Number-of-chunks Total-count

        f.close()
        
    
    def _get_position_frames(self):
        """
        Gets the starting position of each frame (Offsets in bytes), that can be used with seek
        Assumes that the number of chunks is constant
        """
        
        f = open( self.filename)
        f.seek(self._first_line) # start reading after header
        lines_per_frame = self.n_chunks + 1
        offsets = []
        counter = 0
        
        line = True
        while line:
            if not counter % lines_per_frame:
                offsets.append(int(f.tell()))
            line = f.readline()
            counter += 1
        array_offsets = offsets[:-1]  # last is EOF
        self.offsets = np.array(array_offsets, dtype = int)
        
        f.close()

    def __iter__(self):
        """Iterate over time_steps

        ``for chunk in timestep``

            iterate of the chunks
        """
        for i in range(self.n_frames):
            yield self[i]

    def __len__(self):
            return self.n_frames

    def _apply_limits(self, frame):
        if frame < 0:
            frame += len(self)
        if frame < 0 or frame >= len(self):
            raise IndexError("Index {} exceeds length of trajectory ({})."
                             "".format(frame, len(self)))
        return frame


    def _read_frame(self, frame):
        print ("Reading frame %s/%s"%(frame, self.n_frames))
        f = open(self.filename)

        byte_pos = self.offsets[frame:frame+2]
        f.seek(byte_pos[0])
        step, n_chunks, _ = f.readline().split() # Timestep Number-of-chunks Total-count
        start = f.tell()
        stuff = f.read(byte_pos[1] - start)
        data = pd.read_csv(StringIO(stuff),sep=" ",header=None).dropna(axis=1,how='all')
        data.columns = self.header
        self.frame = timestep(frame, n_chunks, data)
        # Example implementation in the DCDReader:
        # self._jump_to_frame(frame)
        # ts = self.ts
        # ts.frame = self._read_next_frame(ts._x, ts._y, ts._z,
        #                                  ts._unitcell, 1)
        return self.frame


    def check_slice_indices(self, start, stop, step):
        """Check frame indices are valid and clip to fit trajectory.

        The usage follows standard Python conventions for :func:`range` but see
        the warning below.

        Parameters
        ----------
        start : int or None
          Starting frame index (inclusive). ``None`` corresponds to the default
          of 0, i.e., the initial frame.
        stop : int or None
          Last frame index (exclusive). ``None`` corresponds to the default
          of n_frames, i.e., it includes the last frame of the trajectory.
        step : int or None
          step size of the slice, ``None`` corresponds to the default of 1, i.e,
          include every frame in the range `start`, `stop`.

        Returns
        -------
        start, stop, step : tuple (int, int, int)
          Integers representing the slice

        Warning
        -------
        The returned values `start`, `stop` and `step` give the expected result
        when passed in :func:`range` but gives unexpected behavior when passed
        in a :class:`slice` when ``stop=None`` and ``step=-1``

        This can be a problem for downstream processing of the output from this
        method. For example, slicing of trajectories is implemented by passing
        the values returned by :meth:`check_slice_indices` to :func:`range` ::

          range(start, stop, step)

        and using them as the indices to randomly seek to. On the other hand,
        in :class:`MDAnalysis.analysis.base.AnalysisBase` the values returned
        by :meth:`check_slice_indices` are used to splice the trajectory by
        creating a :class:`slice` instance ::

          slice(start, stop, step)

        This creates a discrepancy because these two lines are not equivalent::

            range(10, -1, -1)             # [10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0]
            range(10)[slice(10, -1, -1)]  # []

        """

        slice_dict = {'start': start, 'stop': stop, 'step': step}
        for varname, var in slice_dict.items():
            if isinstance(var, numbers.Integral):
                slice_dict[varname] = int(var)
            elif (var is None):
                pass
            else:
                raise TypeError("{0} is not an integer".format(varname))

        start = slice_dict['start']
        stop = slice_dict['stop']
        step = slice_dict['step']

        if step == 0:
            raise ValueError("Step size is zero")

        nframes = len(self)
        step = step or 1

        if start is None:
            start = 0 if step > 0 else nframes - 1
        elif start < 0:
            start += nframes
        if start < 0:
            start = 0

        if step < 0 and start >= nframes:
            start = nframes - 1

        if stop is None:
            stop = nframes if step > 0 else -1
        elif stop < 0:
            stop += nframes

        if step > 0 and stop > nframes:
            stop = nframes

        return start, stop, step

    def __getitem__(self, frame):
        """(Copied from MDANALYSIS, base) 
        Return the Timestep corresponding to *frame*.

        If `frame` is a integer then the corresponding frame is
        returned. Negative numbers are counted from the end.

        If frame is a :class:`slice` then an iterator is returned that
        allows iteration over that part of the trajectory.

        Note
        ----
        *frame* is a 0-based frame index.
        """
        if isinstance(frame, int):
            frame = self._apply_limits(frame)
            return self._read_frame(frame)
        elif isinstance(frame, (list, np.ndarray)):
            if len(frame) != 0 and isinstance(frame[0], (bool, np.bool_)):
                # Avoid having list of bools
                frame = np.asarray(frame, dtype=np.bool)
                # Convert bool array to int array
                frame = np.arange(len(self))[frame]
            return FrameIteratorIndices(self, frame)
        elif isinstance(frame, slice):
            start, stop, step = self.check_slice_indices(
                frame.start, frame.stop, frame.step)
            if start == 0 and stop == len(self) and step == 1:
                return FrameIteratorAll(self)
            else:
                return FrameIteratorSliced(self, frame)
        else:
            raise TypeError("Trajectories must be an indexed using an integer,"
                            " slice or list of indices")
    
        
        #TODO use something as in MDANALYSIS "read_frame", or read next frame, 
        #ir iter, such that when I am iterating I do not open and close the file
        
#    def get_all_data(self):
#        
#        self.frames.append(timestep())
#        f = open( self.filename)
#        f.seek(self._first_line) # start reading after header
#        lines_per_frame = self.n_chunks + 1
#        offsets = []
#        counter = 0
#        
#        line = True
#        while line:
#            if not counter % lines_per_frame:
#                offsets.append(int(f.tell()))
#            line = f.readline()
#            counter += 1
#        array_offsets = offsets[:-1]  # last is EOF
#        self._offsets = np.array(array_offsets, dtype = int)
#        
#        f.close()
        
        
# One way of using it
#results = chunk_reader("velocity_all.dat") 
#        
#    
## Creating all the frames
#
#f = open( results.filename)
#
#series = []
#
#byte_pos = results.offsets
#
#for i, start in enumerate(byte_pos[:-1]):
#
#    f.seek(start)
#    step, n_chunks, _ = f.readline().split() # Timestep Number-of-chunks Total-count
#    start = f.tell()
#    stuff = f.read(byte_pos[i+1] - start)
#    data = pd.read_csv(StringIO(stuff),sep=" ",header=None).dropna(axis=1,how='all')
#    data.columns = results.header
#    series.append(timestep(step, n_chunks, data))


