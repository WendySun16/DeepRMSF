#ifndef chimera_PseudoBondGroup_h
#define	chimera_PseudoBondGroup_h
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

#include <vector>
#include <list>
#include <otf/WrapPy2.h>
#include <_chimera/LineType.h>
#include "_molecule_config.h"
#include <_chimera/Model.h>
#include "PseudoBond.h"
#include "chimtypes.h"

namespace molecule {
class PseudoBondMgr;

class MOLECULE_IMEX PseudoBondGroup: public chimera::Model  {
	friend class PseudoBondMgr;
	void	operator=(const PseudoBondGroup &);	// disable
		PseudoBondGroup(const PseudoBondGroup &);	// disable
protected:
		virtual ~PseudoBondGroup();
	std::vector<PseudoBond *>	PseudoBonds_;
public:
	PseudoBond	*newPseudoBond(Atom *a1, Atom *a2);
	void	deletePseudoBond(PseudoBond *element);
	typedef std::vector<PseudoBond *> PseudoBonds;
	inline const PseudoBonds	&pseudoBonds() const;
	PseudoBond	*findPseudoBond(size_t) const;
public:
//	void	moleculePseudoBonds(std::list<PseudoBond *> *pbl, Molecule *m1,
//	  Molecule *m2 = NULL) const;
//	  // return the PseudoBonds of the given Molecule, or between
//	  // the pair of given Molecules
	
	inline Symbol category() const;
	  // the type of PseudoBonds stored in this group
	
	void deleteAll();
	  // clear the PseudoBondGroup
private:
	Symbol _category;
	void _recategorize(Symbol newCategory);
	struct Reason: public NotifierReason {
                Reason(const char *r): NotifierReason(r) {}
        };
	static Reason		CATEGORY_CHANGED;
public:
	inline bool	showStubBonds() const;
	void		setShowStubBonds(bool);
	inline float	lineWidth() const;
	void		setLineWidth(float);
	inline float	stickScale() const;
	void		setStickScale(float);
	inline LineType	lineType() const;
	void		setLineType(LineType linetype);
	inline void	wireStipple(/*OUT*/ int *factor, /*OUT*/ int *pattern) const;
	void		setWireStipple(int factor, int pattern);
	void		lineTypeWireStipple(int *factor, int *pattern) const;

	// virtual functions from Model
	virtual bool	computeBounds(/*OUT*/ Sphere *s, /*OUT*/ BBox *box) const;
	void		destroy();
	// virtual functions from Selectable
	Selectable::Selectables oslChildren() const;
	Selectable::Selectables oslParents() const;
	Molecule	*allOneMolecule() const;

	virtual PyObject* wpyNew() const;
protected:
	// virtual functions from Model
	LensModel	*newLensModel(Lens *);
	bool		showStubBonds_;
	float		lineWidth_;
	float		stickScale_;
	LineType	lineType_;
	int		stipple[2];
	static TrackChanges::Changes *const
			changes;
	typedef std::set<Molecule *> MolSet;
	mutable MolSet	molecules;
	class MoleculeNotifier: public Notifier {
	public:
                void update(const void *tag, void *style,
                                        const NotifierReason &reason) const;
	}		moleculeNotifier;
	void 		updateMolecule(Molecule *s, const NotifierReason &);
	friend class MoleculeNotifier;
	class OpenModelsNotifier: public Notifier {
	public:
                void update(const void *tag, void *style,
                                        const NotifierReason &reason) const;
	}		openModelsNotifier;
	void 		updateOpenModels(const NotifierReason &);
	friend class OpenModelsNotifier;
	friend class PBGLensModel;
	void		updateModelNotifications() const;
public:
	void		clearModelNotifications();
	virtual void	trackReason(const NotifierReason &reason) const;
	virtual void	wpyAssociate(PyObject* o) const;
protected:
	PseudoBondGroup(PseudoBondMgr *, Symbol category);
};

} // namespace molecule

#endif
