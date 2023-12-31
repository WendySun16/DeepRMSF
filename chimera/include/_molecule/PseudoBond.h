#ifndef chimera_PseudoBond_h
#define	chimera_PseudoBond_h
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

#include <otf/WrapPy2.h>
#include "_molecule_config.h"
#include "Atom.h"
#include "Bond.h"
#include "chimtypes.h"

namespace molecule {
class PseudoBondGroup;

class MOLECULE_IMEX PseudoBond: public Selectable  {
	friend class Atom;
	friend class PseudoBondGroup;
	void	operator=(const PseudoBond &);	// disable
		PseudoBond(const PseudoBond &);	// disable
		virtual ~PseudoBond();
	Array<Atom *, 2>	Atoms_;
	PseudoBondGroup	*PseudoBondGroup_;
public:
	typedef Array<Atom *, 2> Atoms;
	inline const Atoms	&atoms() const;
	Atom	*findAtom(size_t) const;
	PseudoBondGroup	*pseudoBondGroup() const;
public:
	inline Atom	*otherAtom(const Atom *a) const;
	inline bool	contains(const Atom *) const;
private:
	void		basicInitWithAtoms(Atom *, Atom *);
public:
        inline Symbol	category(void) const;
	inline Real		length() const;
	inline Real		sqlength() const;
public:
	typedef Bond::DrawMode DrawMode;
	inline const Color	*color() const;
	void		setColor(/*NULL_OK*/ const Color *);
	inline DrawMode	drawMode() const;
	void		setDrawMode(DrawMode);
	typedef Bond::DisplayMode DisplayMode;
	inline DisplayMode	display() const;
	void		setDisplay(DisplayMode);
	bool		shown() const;
	inline float		radius() const;
	void		setRadius(float);
	inline bool		halfbond() const;
	void		setHalfbond(bool);
	inline const std::string &
			label() const;
	void		setLabel(const std::string &s);
	inline const Vector	&labelOffset() const;
	void		setLabelOffset(const Vector &offset);
	inline const Color	*labelColor() const;
	void		setLabelColor(/*NULL_OK*/ const Color *);
	Point		labelCoord() const;
	Vector		currentLabelOffset() const;
	std::string	oslIdent(Selector start = SelDefault,
					Selector end = SelDefault) const;
	Selectable::Selectables oslChildren() const;
	Selectable::Selectables oslParents() const;
	bool		oslTestAbbr(OSLAbbreviation *a) const;
	inline Selector	oslLevel() const;
	static const Selector	selLevel = SelBond;

	void		reuse(Atom *a0, Atom *A1);

	virtual PyObject* wpyNew() const;
private:
	const Color	*color_;
	DrawMode	drawMode_;
	DisplayMode	display_;
	bool		halfbond_;
	float		radius_;
	std::string	label_;
	Vector		labelOffset_;
	const Color	*labelColor_;
	static TrackChanges::Changes *const
			changes;
	inline void		setMajorChange();
public:
	virtual void	wpyAssociate(PyObject* o) const;

	void		trackReason(const NotifierReason &reason) const;
	static Bond::Reason PSEUDOBOND_REUSED;
private:
	PseudoBond(PseudoBondGroup *, Atom *a1, Atom *a2);
};

} // namespace molecule

#endif
