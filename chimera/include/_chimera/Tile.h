#ifndef Chimera_Tile_h
# define Chimera_Tile_h

# if defined(_MSC_VER) && (_MSC_VER >= 1020)
#  pragma once
# endif

# include "_chimera_config.h"
# include "CameraMode.h"

# ifndef WrapPy

extern "C" {
typedef struct _TRctx TRcontext;
struct Togl;
}

namespace chimera {

class CHIMERA_IMEX Tile: public CameraMode {
	CameraMode	*cm;
	unsigned int	supersample;
	unsigned char	*image;
	unsigned char	*tileImage;
	TRcontext	*trc;
	unsigned int	components;
	Togl		*togl;
	int		start_view;	// starting view in CameraMode
	int		num_views;	// number of views in CameraMode
	int		num_tiles;	// total number of tiles
	double		orthoValues_[6];
public:
	Tile(const char *cameraMode, int imageWidth, int imageHeight,
			int originalWidth, int originalHeight,
			int mode, unsigned char *image, int tileWidth,
			int tileHeight, Togl *t,
			ViewLen viewInfo, unsigned int supersample = 1);
	virtual		~Tile();
	int		numViews() const;
	const CameraView *
			view(int view) const;
	void		computeViews(const Camera &camera);
	bool		setup(const Viewer *viewer, int view) const;
	bool		setupView(const Viewer *viewer, int view, bool ortho) const;
	void		ortho(double left, double right, double bottom,
					double top, double near, double far);
	void		orthoValues(/*OUT*/ double *left, /*OUT*/ double *right,
				/*OUT*/ double *bottom, /*OUT*/ double *top,
				/*OUT*/ double *near, /*OUT*/ double *far);
	void		scissor(int x, int y, int width, int height);
	void		rasterPos3(float x, float y, float z);
	void		rasterPos3(double x, double y, double z);
	void		offsetViewport(int x, int y);
	void		fullscreenRect();
	unsigned int	supersamples();
	void		size(/*OUT*/ int *width, /*OUT*/ int *height);
	void		imageSize(/*OUT*/ int *width, /*OUT*/ int *height);
};

inline unsigned int
Tile::supersamples()
{
	return supersample;
}

} // namespace chimera

# endif /* WrapPy */

#endif
