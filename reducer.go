package reducer

// #cgo CFLAGS: -O3
// #cgo LDFLAGS: -lm
// #include "ikalman/gps.h"
// #include "ikalman/gps.c"
// #include "ikalman/kalman.h"
// #include "ikalman/kalman.c"
// #include "ikalman/matrix.h"
// #include "ikalman/matrix.c"
import "C"

import (
	"math"

	"github.com/paulmach/go.geo"
	"github.com/paulmach/go.geo/reducers"
)

type Reducer struct {
	Smoothing int
	Noise     float64
	Threshold float64
}

func New() *Reducer {
	return &Reducer{
		Smoothing: 10,
		Noise:     200.0,
		Threshold: 1.0,
	}
}

// GeoReduce pulls from a large bag of tricks to remove some of the noise
// that may be found in some slides.
// TODO: Seriously, please improve/simplify this!
func (r *Reducer) GeoReduce(path *geo.Path) *geo.Path {
	geoDistance := path.GeoDistance()
	scaleFactor := geo.MercatorScaleFactor(path.Bound().Center().Lat())
	threshold := r.Threshold * scaleFactor

	path.Transform(geo.Mercator.Project)

	// resample at 1 meter
	path.Resample(int(math.Ceil(geoDistance)))

	// smooth
	smoothed := smooth(path, r.Smoothing)

	// kalman
	kf := C.alloc_filter_velocity2d(C.double(r.Noise))
	for i := 0; i < smoothed.Length()-1; i++ {
		C.update_velocity2d(kf, C.double(smoothed.GetAt(i).Lat()), C.double(smoothed.GetAt(i).Lng()), 1.0)

		var cLat, cLng C.double
		C.get_lat_long(kf, &cLat, &cLng)
		smoothed.SetAt(i, geo.NewPoint(float64(cLng), float64(cLat)))
	}

	smoothed = reducers.VisvalingamThreshold(smoothed, 6*threshold*threshold)
	smoothed = reducers.DouglasPeucker(smoothed, threshold)

	smoothed.Transform(geo.Mercator.Inverse)
	return smoothed
}

// smooth computes a moving average of plus/minus `smooth` points around each point.
func smooth(path *geo.Path, smooth int) *geo.Path {
	smoothed := geo.NewPath()
	smoothed.Push(path.GetAt(0))

	for i := 1; i < path.Length()-1; i++ {
		p := geo.NewPoint(0, 0)
		count := 0

		for j := i - smooth; j <= i+smooth; j++ {
			if j < 0 || j >= path.Length() {
				p.Add(path.GetAt(i))
			} else {
				p.Add(path.GetAt(j))
			}
			count++
		}

		p.Scale(1.0 / float64(count))
		smoothed.Push(p)
	}

	smoothed.Push(path.GetAt(path.Length() - 1))
	return smoothed
}
