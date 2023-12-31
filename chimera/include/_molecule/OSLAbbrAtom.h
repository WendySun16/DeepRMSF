#ifndef Chimera_OSLAbbrAtom_h
# define Chimera_OSLAbbrAtom_h

# if defined(_MSC_VER) && (_MSC_VER >= 1020)
#  pragma once
# endif

# ifndef WrapPy

# include <vector>
# include <string>
# include <pcrecpp.h>
# include <_chimera/Selectable.h>
# include "chimtypes.h"
# include "_molecule_config.h"

namespace molecule {

class Atom;

class MOLECULE_IMEX OSLAbbrAtom : public OSLAbbrTest {

	class Item {
		pcrecpp::RE	*name_;
		std::string	altLoc_[2];
		bool		anyAltLoc_;
		bool		noAltLoc_;
	public:
				Item(const std::string &left,
					const std::string &startRight,
					const std::string &endRight,
					bool hasDot);
				~Item();
		bool		test(const Atom *a) const;
	};
private:
	typedef std::vector<Item *>	Items;
	Items				names_;
public:
			OSLAbbrAtom(const OSLAbbreviation &a);
	virtual		~OSLAbbrAtom();
	virtual void	add(const std::string &left,
				const std::string &right,
				bool hasDot);
	bool		test(const Atom *a) const;
public:
	static std::string	ident(const Atom *a);
};

} // namespace molecule

# endif /* WrapPy */

#endif
