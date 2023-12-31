#ifndef Chimera_LensModel_h
# define Chimera_LensModel_h

# if defined(_MSC_VER) && (_MSC_VER >= 1020)
#  pragma once
# endif

# ifndef WrapPy

# include <set>
# include <utility>

# include "Notifier.h"
# include "LensViewer.h"
# include "Geom3d.h"
# include "Plane.h"

namespace chimera {

// TODO: used hashed set
typedef std::set<ColorGroup *> ColorGroupSet;
typedef std::set<const Color *> ColorSet;

class Model;
class X3DScene;

class CHIMERA_IMEX LensModel: public NotifierList {
	// Per-lens, per-model graphical cache
	virtual void	draw(const LensViewer *viewer, LensViewer::DrawPass pass) const = 0;
	virtual void	drawPick(const LensViewer *viewer) const = 0;
	virtual void	drawPickLabels(const LensViewer *viewer) const;
public:
	void		xformDraw(const LensViewer *viewer, LensViewer::DrawPass pass) const;
	void		xformDrawPick(const LensViewer *viewer) const;
	void		xformDrawPickLabels(const LensViewer *viewer) const;
	void		drawBound(bool sphere = false) const;
	virtual bool	drawLast() const { return false; }
	virtual bool	validXform() const;
	virtual void	invalidateCache() = 0;
	virtual void	invalidateSelection() = 0;
	virtual void	x3dNeeds(/*INOUT*/ X3DScene *scene) const;
	virtual void	x3dWrite(std::ostream &out, unsigned indent,
					/*INOUT*/ X3DScene *scene) const;
	void		destroy();

	// let subclasses know if we're being clipped
	enum Clipped { Hither = 0x1, Cut = 0x2, CutBack = 0x4, Dirty = 0x8 };
	int		isClipped() const;
	void		invalidateClippingPlanes();

	virtual Model	*model() const = 0;

	struct CHIMERA_IMEX Reason: public NotifierReason {
		Reason(const char *r): NotifierReason(r) {}
	};
	static Reason	DIRTY;

	virtual		~LensModel() {}
protected:
	LensModel(Lens *lens);	// should always be on the heap
	Lens		*lens;
	mutable int	clipped;
	mutable Plane	hither, cut, cutback;
};

CHIMERA_IMEX extern void unmonitorColors(ColorSet &colorsUsed,
					 ColorGroupSet &cgsUsed,
					 const void *tag);

inline void
LensModel::destroy()
{
	delete this;
}

inline int
LensModel::isClipped() const
{
	return clipped;
}

} // namespace chimera

# endif /* WrapPy */

#endif
