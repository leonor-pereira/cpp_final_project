// C++ final project header file - PHYS30762

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <cmath>

const double electron_charge = -1.6e-19; // coulombs
const double speed_of_light = 2.99792458e8; // m/s

// mass values
const float up_mass = 2.0;
const float down_mass = 5.0;
const float charm_mass = 1.0*std::pow(10, 3);
const float strange_mass = 95.0;
const float top_mass = 173.0*std::pow(10, 3);
const float bottom_mass = 4.18*std::pow(10, 3);
const float electron_mass = 1.0; // to the nearest mev
const float muon_mass = 106.0; // to the nearest mev
const float tau_mass = 2.0*std::pow(10, 3);
const float higgs_mass = 125.0*std::pow(10, 3);
const float w_mass = 80.0*std::pow(10, 3);
const float z_mass = 91.0*std::pow(10, 3);

// four-momentum class
class four_momentum
{
  private:
    std::shared_ptr<std::vector<float>> four_momentum_vector = std::make_shared<std::vector<float>>(4);

  public:
    // constructors and member functions
    four_momentum() = default;
    four_momentum(float energy, float px, float py, float pz)
    {
      if(energy > 0)
      {
        (*four_momentum_vector).at(0) = energy;
        (*four_momentum_vector).at(1) = px;
        (*four_momentum_vector).at(2) = py;
        (*four_momentum_vector).at(3) = pz;
      }
      else
      {
        std::cout << "The energy value entered is smaller or equal to zero. Setting the energy to 50 GeV. " << std::endl;

        (*four_momentum_vector).at(0) = 50;
        (*four_momentum_vector).at(1) = px;
        (*four_momentum_vector).at(2) = py;
        (*four_momentum_vector).at(3) = pz;
      }
    };
    ~four_momentum()
    {
      std::cout << "Calling destructor..." << std::endl;
    };
    void set_e_px_py_pz(float energy, float px, float py, float pz)
    {
      if(energy > 0)
      {
        (*four_momentum_vector).at(0) = energy;
        (*four_momentum_vector).at(1) = px;
        (*four_momentum_vector).at(2) = py;
        (*four_momentum_vector).at(3) = pz;
      }
      else
      {
        std::cout << "The energy value entered is smaller or equal to zero. Setting the energy to 50 GeV. " << std::endl;
        (*four_momentum_vector).at(0) = 50;
        (*four_momentum_vector).at(1) = px;
        (*four_momentum_vector).at(2) = py;
        (*four_momentum_vector).at(3) = pz;
      }
    };
    void set_e(float energy)
    {
      if(energy > 0)
      {
        (*four_momentum_vector).at(0) = energy;
      }
      else
      {
        std::cout << "The energy value entered is smaller or equal to zero. Setting the energy to 50 GeV. " << std::endl;
        (*four_momentum_vector).at(0) = 50;
      }
    };
    void set_px(float px){(*four_momentum_vector).at(1) = px;};
    void set_py(float py){(*four_momentum_vector).at(2) = py;};
    void set_pz(float pz){(*four_momentum_vector).at(3) = pz;};
    float get_energy() const {return (*four_momentum_vector).at(0);};
    float get_px() const {return (*four_momentum_vector).at(1);};
    float get_py() const {return (*four_momentum_vector).at(2);};
    float get_pz() const {return (*four_momentum_vector).at(3);};

    // copy constructor
    four_momentum(four_momentum &object)
    {
      std::cout << "Calling copy constructor... " << std::endl;

      std::shared_ptr<std::vector<float>> four_momentum_copy = std::make_shared<std::vector<float>>(4);
      set_e_px_py_pz(object.get_energy(), object.get_px(), object.get_py(), object.get_pz());
    };

    // assignment operator
    four_momentum & operator=(const four_momentum &object)
    {
      std::cout << "Calling copy assignment... " << std::endl;

      if(&object == this) {return *this;}
      else
      {
        std::shared_ptr<std::vector<float>> four_momentum_copy = std::make_shared<std::vector<float>>(4);

        set_e_px_py_pz(object.get_energy(), object.get_px(), object.get_py(), object.get_pz());
      }

      return *this;
    };

    // move constructor
    four_momentum(four_momentum &&object) noexcept
    {
      std::cout << "Calling move constructor... " << std::endl;

      (*four_momentum_vector).at(0) = object.get_energy();
      (*four_momentum_vector).at(1) = object.get_px();
      (*four_momentum_vector).at(2) = object.get_py();
      (*four_momentum_vector).at(3) = object.get_pz();

      object.four_momentum_vector = nullptr;
    };

    // move operator
    four_momentum & operator=(four_momentum &&object)
    {
      std::cout << "Calling move assignment... " << std::endl;

      std::swap(four_momentum_vector, object.four_momentum_vector);

      object.four_momentum_vector = nullptr;

      return *this;
    };

    float invariant_mass_calculation()
    {
      float invariant_mass_value_squared, invariant_mass_value;
      invariant_mass_value_squared = std::pow((*four_momentum_vector).at(0),2) - std::pow((*four_momentum_vector).at(1), 2) - std::pow((*four_momentum_vector).at(2), 2) - std::pow((*four_momentum_vector).at(3), 2);
      invariant_mass_value = std::pow(invariant_mass_value_squared, 0.5);

      return invariant_mass_value;
    };
};

// abstract base particle class 
class particle
{
  protected:
    float mass;
    float charge;
    float spin;
    std::string flavour;
    float baryon_number;
    float lepton_number;
    four_momentum four_momentum_vector;
    bool antimatter_status;

  public:
    particle() = default;
    virtual ~particle() { };
    void set_mass(float mass_input)
    {
      if(mass_input >= 0)
      {
        mass = mass_input;
      }
    };
    virtual void set_mass() = 0;
    virtual void set_charge(float charge_input) = 0;
    virtual void set_flavour(std::string flavour_input) = 0;
    void set_four_momenta(four_momentum particle_four_vector) {four_momentum_vector = particle_four_vector;};
    float get_mass() {return mass;};
    float get_charge() {return charge;};
    float get_spin() {return spin;};
    std::string get_flavour() {return flavour;};
    float get_baryon_number() {return baryon_number;};
    float get_lepton_number() {return lepton_number;};
    four_momentum get_four_momenta() {return four_momentum_vector;};
    bool get_antimatter_status() {return antimatter_status;};

    friend std::vector<float> operator+(four_momentum particle_1, four_momentum particle_2)
    {
      std::vector<float> summed_four_vec(4);

      summed_four_vec.at(0) = particle_1.get_energy() + particle_2.get_energy();
      summed_four_vec.at(3) = particle_1.get_pz() + particle_2.get_pz();
      summed_four_vec.at(1) = particle_1.get_px() + particle_2.get_px();
      summed_four_vec.at(2) = particle_1.get_py() + particle_2.get_py();
      
      return summed_four_vec;
    };

    friend std::vector<float> operator-(four_momentum particle_1, four_momentum particle_2)
    {
      std::vector<float> summed_four_vec(4);

      summed_four_vec.at(0) = particle_1.get_energy() - particle_2.get_energy();
      summed_four_vec.at(1) = particle_1.get_px() - particle_2.get_px();
      summed_four_vec.at(2) = particle_1.get_py() - particle_2.get_py();
      summed_four_vec.at(3) = particle_1.get_pz() - particle_2.get_pz();
      
      return summed_four_vec;
    };

    friend float dot_product(four_momentum particle_1, four_momentum particle_2){
      float dot_product_value;

      dot_product_value = (particle_1.get_energy()*particle_2.get_energy()) - (((particle_1.get_px())*(particle_2.get_px())) + ((particle_1.get_py())*(particle_2.get_py())) + ((particle_1.get_pz())*(particle_2.get_pz())));
      return dot_product_value;
    };

    virtual void particle_printing_function() = 0;       //   p00pp00
};

// derived classes
// fermion class
class fermion: public particle
{
  public:
    fermion() = default;
    ~fermion(){ };
    // member functions
    void particle_printing_function() override 
    {
      std::cout << "Printing four-momenta elements. Energy:" << four_momentum_vector.get_energy() << "; px: " << four_momentum_vector.get_px() << "; py: " << four_momentum_vector.get_py() << "; pz: " << four_momentum_vector.get_pz() << "." << std::endl;
      std::cout << "Charge: " << get_charge() << "; spin: " << get_spin() << ". " << std::endl;
    };
};

// classes derived from fermion
class quark: public fermion
{
  private:
    std::string colour;
  public:
    quark() = default;
    quark(std::string flavour_input)
    {
      //spin = 0.5;
      if(flavour_input == "up" || flavour_input == "charm" || flavour_input == "top")
      {
        flavour = flavour_input;
        charge = 2.0/3.0;
        baryon_number = 1.0/3.0;
        lepton_number= 0;
        antimatter_status = 0;
        spin = 1.0/2.0;
      }
      else if(flavour_input == "down" || flavour_input == "strange" || flavour_input == "bottom")
      {
        flavour = flavour_input;
        charge = -1.0/3.0;
        baryon_number = 1.0/3.0;
        lepton_number= 0;
        antimatter_status = 0;
        spin = 1.0/2.0;
      }
      else if(flavour_input == "antiup" || flavour_input == "anticharm" || flavour_input == "antitop")
      {
        flavour = flavour_input;
        charge = -2.0/3.0;
        baryon_number = -1.0/3.0;
        lepton_number= 0;
        antimatter_status = 1;
        spin = 1.0/2.0;
      }
      else if(flavour_input == "antidown" || flavour_input == "antistrange" || flavour_input == "antibottom")
      {
        flavour = flavour_input;
        charge = 1.0/3.0;
        baryon_number = -1.0/3.0;
        lepton_number= 0;
        antimatter_status = 1;
        spin = 1.0/2.0;
      }
      else
      {
        std::cout << "Flavour entered is invalid. Setting flavour to up. " << std::endl;
        flavour = "up";
        charge = 2.0/3.0;
        baryon_number = 1.0/3.0;
        lepton_number= 0;
        antimatter_status = 0;
        spin = 1.0/2.0;
      }
    };
    ~quark(){ };

    // member functions
    void set_flavour(std::string flavour_input) override
    {
      if(flavour_input == "up" || flavour_input == "charm" || flavour_input == "top")
      {
        flavour = flavour_input;
        charge = 2.0/3.0;
        baryon_number = 1.0/3.0;
        lepton_number= 0;
        antimatter_status = 0;
        spin = 1.0/2.0;
      }
      else if(flavour_input == "down" || flavour_input == "strange" || flavour_input == "bottom")
      {
        flavour = flavour_input;
        charge = -1.0/3.0;
        baryon_number = 1.0/3.0;
        lepton_number= 0;
        antimatter_status = 0;
        spin = 1.0/2.0;
      }
      else if(flavour_input == "antiup" || flavour_input == "anticharm" || flavour_input == "antitop")
      {
        flavour = flavour_input;
        charge = -2.0/3.0;
        baryon_number = -1.0/3.0;
        lepton_number= 0;
        antimatter_status = 1;
        spin = 1.0/2.0;
      }
      else if(flavour_input == "antidown" || flavour_input == "antistrange" || flavour_input == "antibottom")
      {
        flavour = flavour_input;
        charge = 1.0/3.0;
        baryon_number = -1.0/3.0;
        lepton_number= 0;
        antimatter_status = 1;
        spin = 1.0/2.0;
      }
      else
      {
        std::cout << "Flavour entered is invalid. Setting flavour to up. " << std::endl;
        flavour = "up";
        charge = 2.0/3.0;
        baryon_number = 1.0/3.0;
        lepton_number= 0;
        antimatter_status = 0;
        spin = 1.0/2.0;
      }
    };

    void set_charge(float charge_input) override
    {
      spin = 1.0/2.0;
      if(charge_input == 2.0/3.0 || charge_input == -1.0/3.0)
      {
        baryon_number = 1.0/3.0;
        lepton_number= 0;
        antimatter_status = 0;
      }
      else if(charge_input == -2.0/3.0 || charge_input == 1.0/3.0)
      {
        baryon_number = -1.0/3.0;
        lepton_number= 0;
        antimatter_status = 1;
      }
      else
      {
        std::cout << "The charge value inputted is invalid. Setting charge to +2/3." << std::endl;
        charge == 2.0/3.0;
        baryon_number = 1.0/3.0;
        lepton_number= 0;
        antimatter_status = 0;
      }
    };

    void set_colour(std::string colour_input)
    {
      spin = 1.0/2.0;

      if((colour_input == "red" || colour_input == "blue" || colour_input == "green") && antimatter_status == 0)
      {
        colour = colour_input;
      }
      else if((colour_input == "antired" || colour_input == "antiblue" || colour_input == "antigreen") && antimatter_status == 1)
      {
        colour = colour_input;
      }
      else
      {
        std::cout << "Colour inputted is invalid. Setting colour to red. " << std::endl;
        colour = "red";
      }
    };

    void set_mass() override
    {
      float invariant_mass_value = four_momentum_vector.invariant_mass_calculation();

      if((flavour == "up" || flavour == "antiup") && (invariant_mass_value == up_mass)) {mass = invariant_mass_value;}
      else if((flavour == "down" || flavour == "antidown") && (invariant_mass_value == down_mass)) {mass = invariant_mass_value;}
      else if((flavour == "charm" || flavour == "anticharm") && (invariant_mass_value == charm_mass)) {mass = invariant_mass_value;}
      else if((flavour == "strange" || flavour == "antistrange") && (invariant_mass_value == strange_mass)) {mass = invariant_mass_value;}
      else if((flavour == "top" || flavour == "antitop") && (invariant_mass_value == top_mass)) {mass = invariant_mass_value;}
      else if((flavour == "bottom" || flavour == "antibottom") && (invariant_mass_value == bottom_mass)) {mass = invariant_mass_value;}
      else
      {
        float new_invariant_mass;

        if(flavour == "up" || flavour == "antiup")
        {
          std::cout << "The invariant mass calculated from this particle's four-momentum does not match the expected mass value. Scaling the four-momentum. " << std::endl;
          four_momentum_vector.set_e_px_py_pz(4, 2, 2, 2);
          new_invariant_mass = four_momentum_vector.invariant_mass_calculation();
          mass = new_invariant_mass;
        }
        else if(flavour == "down" || flavour == "antidown")
        {
          std::cout << "The invariant mass calculated from this particle's four-momentum does not match the expected mass value. Scaling the four-momentum. " << std::endl;
          four_momentum_vector.set_e_px_py_pz(10, 5, 5, 5);
          new_invariant_mass = four_momentum_vector.invariant_mass_calculation();
          mass = new_invariant_mass;
        }
        else if(flavour == "charm" || flavour == "anticharm")
        {
          std::cout << "The invariant mass calculated from this particle's four-momentum does not match the expected mass value. Scaling the four-momentum. " << std::endl;
          four_momentum_vector.set_e_px_py_pz(1405, 720, 377, 560);
          new_invariant_mass = four_momentum_vector.invariant_mass_calculation();
          mass = new_invariant_mass;
        }
        else if(flavour == "strange" || flavour == "antistrange")
        {
          std::cout << "The invariant mass calculated from this particle's four-momentum does not match the expected mass value. Scaling the four-momentum. " << std::endl;
          four_momentum_vector.set_e_px_py_pz(108, 37, 11, 34);
          new_invariant_mass = four_momentum_vector.invariant_mass_calculation();
          mass = new_invariant_mass;
        }
        else if(flavour == "top" || flavour == "antitop")
        {
          std::cout << "The invariant mass calculated from this particle's four-momentum does not match the expected mass value. Scaling the four-momentum. " << std::endl;
          four_momentum_vector.set_e_px_py_pz(223158, 63500, 94600, 83000);
          new_invariant_mass = four_momentum_vector.invariant_mass_calculation();
          mass = new_invariant_mass;
        }
        else if(flavour == "bottom" || flavour == "antibottom")
        {
          std::cout << "The invariant mass calculated from this particle's four-momentum does not match the expected mass value. Scaling the four-momentum. " << std::endl;
          four_momentum_vector.set_e_px_py_pz(4941, 1401, 1900, 1170);
          new_invariant_mass = four_momentum_vector.invariant_mass_calculation();
          mass = new_invariant_mass;
        }
        else
        {
          std::cout << "No valid flavour has been set to this particle. Setting flavour to up and setting an appropriate four-momenta. " << std::endl;
          flavour = "up";
          four_momentum_vector.set_e_px_py_pz(4, 2, 2, 2);
          new_invariant_mass = four_momentum_vector.invariant_mass_calculation();
          mass = new_invariant_mass;
        }
      }
    };

    std::string get_colour() {return colour;}

    void particle_printing_function() override 
    {
      std::cout << "Printing quark information." << std::endl;
      fermion::particle_printing_function();
      std::cout << "Quark flavour: " << get_flavour() << "; mass: " << get_mass() << "; colour:" << get_colour() << "; baryon number: " << get_baryon_number() << ". \n" << std::endl;
    };
};

class lepton: public fermion
{
  public:
    lepton() = default;
    ~lepton(){ };
    // member functions
    void particle_printing_function() override 
    {
      fermion::particle_printing_function();
      std::cout << "Lepton number: " << get_lepton_number() <<  std::endl;
    };
};

// classes derived from lepton
class electron: public lepton
{
  private:
    std::shared_ptr<std::vector<float>> energy_deposited = std::make_shared<std::vector<float>>(4);
  public:
    electron() = default;
    ~electron(){ };

    void set_energy_deposited(four_momentum four_vec, float em_1_deposited, float em_2_deposited, float had_1_deposited, float had_2_deposited)
    {
      float energy_sum;
      energy_sum = em_1_deposited + em_2_deposited + had_1_deposited + had_2_deposited;

      if (energy_sum == four_vec.get_energy())
      {
        (*energy_deposited).at(0) = em_1_deposited;
        (*energy_deposited).at(1) = em_2_deposited;
        (*energy_deposited).at(2) = had_1_deposited;
        (*energy_deposited).at(3) = had_2_deposited;
      }
      else
      {
        std::cout << "The four vector entered and the deposited energies do not match. Correcting values. " << std::endl;
        float correcting_factor = (four_vec.get_energy())/energy_sum;
        (*energy_deposited).at(0) = em_1_deposited*correcting_factor;
        (*energy_deposited).at(1) = em_2_deposited*correcting_factor;
        (*energy_deposited).at(2) = had_1_deposited*correcting_factor;
        (*energy_deposited).at(3) = had_2_deposited*correcting_factor;
      }
    };

    void set_flavour(std::string flavour_input) override
    {
      if(flavour_input == "electron")
      {
        flavour = flavour_input;
        charge = -1.0;
        baryon_number = 0.0;
        lepton_number = 1;
        antimatter_status = 0;
        spin = 1.0/2.0;
      }
      else if(flavour_input == "antielectron")
      {
        flavour = flavour_input;
        charge = 1.0;
        baryon_number = 0.0;
        lepton_number = -1;
        antimatter_status = 1;
        spin = 1.0/2.0;
      }
      else
      {
        std::cout << "The flavour inputted is invalid. Setting flavour to electron. " << std::endl;
        flavour = "electron";
        charge = -1.0;
        baryon_number = 0.0;
        lepton_number = 1;
        antimatter_status = 0;
        spin = 1.0/2.0;
      }
    };

    void set_charge(float charge_input) override
    {
      if (charge_input==-1)
      {
        charge = charge_input;
        flavour = "electron";
        baryon_number = 0.0;
        lepton_number = 1;
        antimatter_status = 0;
        spin = 1.0/2.0;
      }
      else if (charge_input==+1)
      {
        charge = charge_input;
        flavour = "antielectron";
        baryon_number = 0.0;
        lepton_number = -1;
        antimatter_status = 1;
        spin = 1.0/2.0;
      }
      else
      {
        std::cout << "Charge value inputted is invalid. Setting charge to -1. " << std::endl;
        charge = -1.0;
        flavour = "electron";    
        baryon_number = 0.0;
        lepton_number = 1;
        antimatter_status = 0;
        spin = 1.0/2.0;    
      }
    };

    std::shared_ptr<std::vector<float>> get_energy_deposited() {return energy_deposited;};

    void particle_printing_function() override 
    {
      std::cout << "Printing electron information. " << std::endl;
      lepton::particle_printing_function();
      std::cout << "Flavour: " << get_flavour() << ". " << "Energy deposited in calorimeter. Electromagnetic 1st layer component: " << (*energy_deposited).at(0) << "; electromagnetic 2nd layer component: " << (*energy_deposited).at(1) << "; hadronic 1st layer component: " << (*energy_deposited).at(2) << "; hadronic 2nd layer component: " << (*energy_deposited).at(3) << ". \n" <<  std::endl;
    };
};

class muon: public lepton
{
  private:
    bool isolation;
  public:
    muon() = default;
    ~muon(){ };
    void set_isolation(bool isolation_value)
    {
      if (isolation_value == 1)
      {
        isolation = isolation_value;
        std::cout << "The muon is isolated. " << std::endl;
      }
      else if (isolation_value == 0)
      {
        isolation = isolation_value;
        std::cout << "The muon is not isolated. " << std::endl;
      }
      else
      {
        std::cout << "An invalid value was inputted. Setting isolation to true. " << std::endl;
        isolation = 1;
      }
    };

    void set_charge(float charge_input) override
    {
      if (charge_input==-1)
      {
        charge = charge_input;
        flavour = "muon";
        baryon_number = 0.0;
        lepton_number = 1;
        antimatter_status = 0;
        spin = 1.0/2.0;
      }
      else if (charge_input==+1)
      {
        charge = charge_input;
        flavour = "antimuon";
        baryon_number = 0.0;
        lepton_number = -1;
        antimatter_status = 1;
        spin = 1.0/2.0;
      }
      else
      {
        std::cout << "Charge value inputted is invalid. Setting charge to -1. " << std::endl;
        charge = -1;
        flavour = "muon";
        baryon_number = 0.0;
        lepton_number = 1;
        antimatter_status = 0;
        spin = 1.0/2.0;   
      }
    };

    void set_flavour(std::string flavour_input) override
    {
      if (flavour_input == "muon")
      {
        flavour = flavour_input;
        charge = -1;
        flavour = "muon";
        baryon_number = 0.0;
        lepton_number = 1;
        antimatter_status = 0;
        spin = 1.0/2.0;
      }
      else if (flavour_input == "antimuon")
      {
        flavour = flavour_input;
        charge = +1;
        baryon_number = 0.0;
        lepton_number = -1;
        antimatter_status = 1;
        spin = 1.0/2.0;
      }
      else
      {
        std::cout << "Flavour inputted is invalid. Setting flavour to muon. " << std::endl;
        charge = -1;
        flavour = "muon";
        baryon_number = 0.0;
        lepton_number = 1;
        antimatter_status = 0;
        spin = 1.0/2.0;
      }
    };

    bool get_isolation() const {return isolation;};

    void particle_printing_function() override 
    {
      if(isolation == 0)
      {
        std::cout << "Printing muon information. " << std::endl;
        lepton::particle_printing_function();
        std:: cout << "Flavour: " << get_flavour() << "; isolation status: false." << std::endl;
      }
      else if(isolation == 1)
      {
        std::cout << "Printing muon information. " << std::endl;
        lepton::particle_printing_function();
        std:: cout << "Flavour: " << get_flavour() << "; isolation status: true." << std::endl;
      }
      else
      {
        std::cout << "Particle does not have all attributes associated. " << std::endl;
      }
    };
};

class tau: public lepton
{
  public:
    tau() = default;
    ~tau(){ };
    //member functions
    void set_flavour(std::string flavour_input) override
    {
      if (flavour_input == "tau")
      {
        flavour = flavour_input;
        charge = -1;
        flavour = "muon";
        baryon_number = 0.0;
        lepton_number = 1;
        antimatter_status = 0;
        spin = 1.0/2.0;
      }
      else if (flavour_input == "antitau")
      {
        flavour = flavour_input;
        charge = +1;
        baryon_number = 0.0;
        lepton_number = -1;
        antimatter_status = 1;
        spin = 1.0/2.0;
      }
      else
      {
        std::cout << "Flavour inputted is invalid. Setting flavour to tau. " << std::endl;
        charge = -1;
        flavour = "tau";
        baryon_number = 0.0;
        lepton_number = 1;
        antimatter_status = 0;
        spin = 1.0/2.0;
      }
    };

    void set_charge(float charge_input) override
    {
      if (charge_input==-1)
      {
        charge = charge_input;
        flavour = "tau";
        baryon_number = 0.0;
        lepton_number = 1;
        antimatter_status = 0;
        spin = 1.0/2.0;
      }
      else if (charge_input==+1)
      {
        charge = charge_input;
        flavour = "antitau";
        baryon_number = 0.0;
        lepton_number = -1;
        antimatter_status = 1;
        spin = 1.0/2.0;
      }
      else
      {
        std::cout << "Charge value inputted is invalid. Setting charge to -1. " << std::endl;
        charge = -1;
        flavour = "tau";
        baryon_number = 0.0;
        lepton_number = 1;
        antimatter_status = 0;
        spin = 1.0/2.0;   
      }
    };
};

class neutrino: public lepton
{
  private:
    bool hasInteracted;
  public:
    neutrino() = default;
    ~neutrino(){ };
    void set_neutrino_flavour(std::string neutrino_flavour)
    {
      charge = 0;
      if (neutrino_flavour == "electron_neutrino" || neutrino_flavour == "muon_neutrino" || neutrino_flavour == "tau_neutrino")
      {
        flavour = neutrino_flavour;
        baryon_number = 0.0;
        lepton_number = 1;
        antimatter_status = 0;
        spin = 1.0/2.0;
      }
      else if(neutrino_flavour == "antielectron_neutrino" || neutrino_flavour == "antimuon_neutrino" || neutrino_flavour == "antitau_neutrino")
      {
        flavour = neutrino_flavour;
        baryon_number = 0.0;
        lepton_number = -1;
        antimatter_status = 1;
        spin = 1.0/2.0;
      }
      else
      {
        std::cout << "The flavour entered is invalid. Setting neutrino flavour to electron neutrino. " << std::endl;
        flavour = "electron_neutrino";
        baryon_number = 0.0;
        lepton_number = 1;
        antimatter_status = 0;
        spin = 1.0/2.0;
      }
    };

    void set_interaction(bool interaction_value)
    {
      if (interaction_value == 1)
      {
        hasInteracted = interaction_value;
        std::cout << "The neutrino has interacted. " << std::endl;
      }
      else if (interaction_value == 0)
      {
        hasInteracted = interaction_value;
        std::cout << "The neutrino has not interacted. " << std::endl;
      }
      else
      {
        std::cout << "An invalid value was inputted. Setting interaction to true. " << std::endl;
        hasInteracted = 1;
      }
    };

    bool get_interaction() const {return hasInteracted;};

    void particle_printing_function() override 
    {
      if(hasInteracted == 0)
      {
        std::cout << "Printing neutrino information. " << std::endl;
        lepton::particle_printing_function();
        std:: cout << "Flavour: " << get_flavour() << " neutrino; interaction status: false." << std::endl;
      }
      else if(hasInteracted == 1)
      {
        std::cout << "Printing neutrino information. " << std::endl;
        lepton::particle_printing_function();
        std:: cout << "Flavour: " << get_flavour() << " neutrino; interaction status: true." << std::endl;
      }
      else
      {
        std::cout << "Particle does not have all attributes associated. " << std::endl;
      }
    };
};

// scalar bosons class
class scalar_boson: public particle
{
  public:
    scalar_boson() = default;
    ~scalar_boson(){ };
    // member functions
    void particle_printing_function() override 
    {
      std::cout << "Printing four-momenta elements. Energy:" << four_momentum_vector.get_energy() << "; px: " << four_momentum_vector.get_px() << "; py: " << four_momentum_vector.get_py() << "; pz: " << four_momentum_vector.get_pz() << "." << std::endl;
      std::cout << "Charge: " << get_charge() << "; spin: " << get_spin() << ". " << std::endl;
    };
};

// classes derived from scalar bosons
class higgs_boson: public scalar_boson
{
  public:
    higgs_boson() = default;
    ~higgs_boson(){ };
    // member functions
    void particle_printing_function() override 
    {
      std::cout << "Printing Higgs boson information. " << std::endl;
      scalar_boson::particle_printing_function();
      // INCLUDE DECAYS HERE
    }
};

// vector bosons class
class vector_boson: public particle
{
  public:
    vector_boson() = default;
    ~vector_boson(){ };
    // member functions
    void particle_printing_function() override 
    {
      std::cout << "Printing four-momenta elements. Energy:" << four_momentum_vector.get_energy() << "; px: " << four_momentum_vector.get_px() << "; py: " << four_momentum_vector.get_py() << "; pz: " << four_momentum_vector.get_pz() << "." << std::endl;
      std::cout << "Charge: " << get_charge() << "; spin: " << get_spin() << ". " << std::endl;
    };
};

// classes derived from massless bosons
class gluon: public vector_boson
{
  private:
    std::vector<std::string> flavour_combination;
  public:
    gluon() = default;
    ~gluon(){ };
    // member functions
    void set_gluon_colour(std::string flavour1, std::string flavour2)
    {
      flavour = "gluon";
      charge = 0.0;
      spin = 1.0;
      baryon_number = 0;
      lepton_number = 0;

      if((flavour1 == "red" && flavour2 == "antiblue") || (flavour1 == "red" && flavour2 == "antigreen") || (flavour1 == "blue" && flavour2 == "antired") || (flavour1 == "blue" && flavour2 == "antigreen") || (flavour1 == "green" && flavour2 == "antired") || (flavour1 == "green" && flavour2 == "antiblue"))
      {
        flavour_combination.at(0) = flavour1;
        flavour_combination.at(1) = flavour2;
      }
      else
      {
        std::cout << "The colour combination entered is invalid. Setting the gluon's colours red and antiblue. " << std::endl;
        flavour_combination.at(0) = "red";
        flavour_combination.at(1) = "antiblue"; 
      }
    };

    std::vector<std::string> get_gluon_colour() {return flavour_combination;};

    void particle_printing_function() override 
    {
      std::cout << "Printing gluon information. " << std::endl;
      vector_boson::particle_printing_function();
      std::cout << "Colour combination: " << flavour_combination.at(0) << " and " << flavour_combination.at(1) << ". " << std::endl;
    };
};

class photon: public vector_boson
{
  public:
    photon() = default;
    ~photon(){ };
    // member functions
    void particle_printing_function() override 
    {
      std::cout << "Printing photon information. " << std::endl;
      vector_boson::particle_printing_function();
    };
};

// classes derived from massive bosons
class w_boson: public vector_boson
{
  public:
    w_boson() = default;
    ~w_boson(){ };
    // member functions
    void set_flavour(std::string flavour_input)
    {
      if(flavour_input == "w_plus")
      {
        flavour = flavour_input;
        charge = +1.0;
        spin = 1.0;
        baryon_number = 0;
        lepton_number = 0;
      }
      else if(flavour_input == "w_minus")
      {
        flavour = flavour_input;
        charge = -1.0;
        spin = 1.0;
        baryon_number = 0;
        lepton_number = 0;
      }
      else
      {
        std::cout << "Flavour inputted is invalid. Setting flavour to W plus. " << std::endl;
        flavour = "w_plus";
        charge = +1.0;
        spin = 1.0;
        baryon_number = 0;
        lepton_number = 0;
      }
    };

    void set_charge(float charge_input) override
    {
      if(charge_input == +1)
      {
        charge = charge_input;
        flavour = "w_plus";
        spin = 1.0;
        baryon_number = 0;
        lepton_number = 0;
      }
      else if(charge_input == -1)
      {
        charge = charge_input;
        flavour = "w_plus";
        spin = 1.0;
        baryon_number = 0;
        lepton_number = 0;        
      }
    };
};

class z_boson: public vector_boson
{
  public:
    z_boson() = default;
    ~z_boson(){ };
};