import csv
import pickle


class UserPathways:
    def __init__(self,
                 user_pw_file,
                 db_name,
                 base_dir,
                 custom_pw_name,
                 custom_pw_exists):
        self.base_dir = base_dir
        self.custom_pw_name = custom_pw_name

        if custom_pw_exists == 'y':
            self.pickled_file = '{}/databases/user_pathways/{}.pkl'.format(
                self.base_dir,
                self.custom_pw_name)
            self.combined_user_pathways_and_db = pickle.load(open(self.pickled_file, 'rb'))
        else:
            self.user_pws = self.read_in_user_pw_file(user_pw_file)
            self.db_name = db_name
            self.db = pickle.load(open('{}/databases/{}.pkl'.format(self.base_dir, self.db_name), 'rb'))
            self.combined_user_pathways_and_db = self.combine_user_pathways_and_db()
            self.pickled_file = self.pickle_custom_pw()

    @staticmethod
    def read_in_user_pw_file(user_pw_file):
        """
        user_pws = {
            all: set(),
            dict: {
                pw_name: {
                    db: str,
                    genes: set,
                    metabolites: set,
                    drugs: set
                }
            }
        }
        """
        user_pws = {
            'all': set(),
            'dict': {}
        }

        with open(user_pw_file, 'r') as tsvin:
            tsvin = csv.reader(tsvin, delimiter=",")

            for row in tsvin:
                pw_name = row[0]
                db = row[1]
                genes = set(row[2:])
                metabolites = set()
                drugs = set()
                user_pws['dict'][pw_name] = {
                    'db': db,
                    'genes': genes,
                    'metabolites': metabolites,
                    'drugs': drugs
                }
                user_pws['all'] = user_pws['all'].union(genes)

        return user_pws

    def combine_user_pathways_and_db(self):
        combined_all = self.db['all'].union(self.user_pws['all'])
        combined_dict = {**self.db['dict'], **self.user_pws['dict']}

        return {
            'all': combined_all,
            'dict': combined_dict
        }

    def pickle_custom_pw(self):
        pickled_file = '{}/databases/user_pathways/{}.pkl'.format(self.base_dir,
                                                                  self.custom_pw_name)
        pickle.dump(self.combined_user_pathways_and_db,
                    open(pickled_file, 'wb'))

        return pickled_file
